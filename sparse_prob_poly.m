function [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_poly(A,b,x_k,U,param)
% convert the polynomial problem with potytope constraints into sparsePOP format

% -------------------- %
% To check the objPoly.coef, run the following after running the
% *_example.m file with the symbolic section uncommented.
% It symbolically finds the expression for \delta \xi.
%{
syms x1 x2 u x_k1 x_k2 gamma
xi = @(x) dhx(x)*fx(x) + dhx(x)*gx(x)*u + gamma*hx(x);
delta_xi = xi([x1,x2])-xi([x_k1,x_k2]);
[c,t] = coeffs(delta_xi, [x1 x2 u]);
%}
% c shows the coefs, t shows the corresponding term types
% Make sure to comment the symbolic section in the *_example.m file after 
% this check.
% -------------------- %

% objPoly 

% Note that the objective polynomial is changed to the equivalent one with positive variables of lower 
% bound 0 and free variables
if strcmp(param.example,'jankovic')
    objPoly.typeCone = 1;
    objPoly.dimVar   = 3;
    objPoly.degree   = 4;
    if param.d-param.gamma == 0
        objPoly.noTerms  = 6;
        objPoly.supports = [0,0,0; 3,1,0; 0,1,0; 0,2,0; 0,2,1; 0,0,1];
        objPoly.coef     = [- (param.q*2*x_k(2)*x_k(1)^3+x_k(2)+param.d*x_k(1)+param.gamma*(param.q*x_k(2)^2-x_k(1)+1)) + param.gamma...
            ; param.q*2; 1; param.gamma*param.q; param.q*2; -param.q*2*x_k(2)^2];
    %     objPoly.coef     = [- (param.q*2*x_k(2)*x_k(1)^3+x_k(2)+param.d*x_k(1)+param.gamma*(param.q*x_k(2)^2-x_k(1)+1)) + param.gamma + param.q*2*x_k(2)^2 ...
    %         ; param.q*2; 1; param.gamma*param.q-param.q*2; param.q*2; -param.q*2*x_k(2)^2];
    else
        objPoly.noTerms  = 7;
        objPoly.supports = [0,0,0; 3,1,0; 0,1,0; 0,2,0; 0,2,1; 0,0,1; 1,0,0];
        objPoly.coef     = [- (param.q*2*x_k(2)*x_k(1)^3+x_k(2)+param.d*x_k(1)+param.gamma*(param.q*x_k(2)^2-x_k(1)+1)) + param.gamma...
            ; param.q*2; 1; param.gamma*param.q; param.q*2; -param.q*2*x_k(2)^2; param.d-param.gamma];
    %     objPoly.coef     = [- (param.q*2*x_k(2)*x_k(1)^3+x_k(2)+param.d*x_k(1)+param.gamma*(param.q*x_k(2)^2-x_k(1)+1)) + param.gamma + param.q*2*x_k(2)^2 ...
    %         ; param.q*2; 1; param.gamma*param.q-param.q*2; param.q*2; -param.q*2*x_k(2)^2; param.d-param.gamma];
    end
elseif strcmp(param.example,'spring')
    objPoly.typeCone = 1;
    objPoly.dimVar   = 3;
    objPoly.degree   = 2;
    objPoly.noTerms  = 8;
    objPoly.supports = [2,0,0; 1,1,0; 1,0,1; 1,0,0; 0,2,0; 0,1,0; 0,0,1; 0,0,0];
    objPoly.coef     = [ param.k^2/param.mass - (param.gamma0*param.gamma*param.k)/2; (param.b*param.k)/param.mass - param.gamma*param.k - param.gamma0*param.k ...
        ; -param.k/param.mass; param.gamma0*param.gamma*param.k*param.x1_shift - (param.k^2*param.x1_shift)/param.mass; - param.k; param.gamma0*param.k*param.x1_shift + param.gamma*param.k*param.x1_shift - (param.b*param.k*param.x1_shift)/param.mass ...
        ; (param.k*param.x1_shift)/param.mass - (param.k*(param.x1_shift - x_k(1)))/param.mass ...
        ; param.gamma0*param.gamma*(param.P_max - (param.k*param.x1_shift^2)/2) - param.gamma*(param.gamma0*(param.P_max - (param.k*(param.x1_shift - x_k(1))^2)/2) + param.k*x_k(2)*(param.x1_shift - x_k(1))) + ...
        param.k*x_k(2)*(x_k(2) - param.gamma0*(param.x1_shift - x_k(1))) + param.k*((param.b*x_k(2))/param.mass + (param.k*x_k(1))/param.mass)*(param.x1_shift - x_k(1))];
else
    error([param.example,' is not a valid example']);
end

%%% polytope bounds
% ineqPolySys (Ax<=b)                        
ineqPolySys = cell(size(b,1),1);                        
for i = 1:size(b,1)
    ineqPolySys{i}.typeCone = 1;
    ineqPolySys{i}.dimVar   = 3;
    ineqPolySys{i}.degree   = 1;
    if A(i,1) == 0
        ineqPolySys{i}.noTerms  = 2;
        ineqPolySys{i}.supports = [0,0,0; 0,1,0];
        ineqPolySys{i}.coef     = [b(i); -A(i,2)];
    elseif A(i,2) == 0
        ineqPolySys{i}.noTerms  = 2;
        ineqPolySys{i}.supports = [0,0,0; 1,0,0];
        ineqPolySys{i}.coef     = [b(i); -A(i,1)];
    else
        ineqPolySys{i}.noTerms  = 3;
        ineqPolySys{i}.supports = [0,0,0; 1,0,0; 0,1,0];
        ineqPolySys{i}.coef     = [b(i); -A(i,1); -A(i,2)];
    end
end

% U(1)<=u<=U(2)
ineqPolySys{i+1}.typeCone = 1;
ineqPolySys{i+1}.dimVar   = 3;
ineqPolySys{i+1}.degree   = 1;
ineqPolySys{i+1}.noTerms  = 2;
ineqPolySys{i+1}.supports = [0,0,0; 0,0,1];
ineqPolySys{i+1}.coef     = [-U(1); 1];

ineqPolySys{i+2}.typeCone = 1;
ineqPolySys{i+2}.dimVar   = 3;
ineqPolySys{i+2}.degree   = 1;
ineqPolySys{i+2}.noTerms  = 2;
ineqPolySys{i+2}.supports = [0,0,0; 0,0,1];
ineqPolySys{i+2}.coef     = [U(2); -1];

% lower bounds for variables x1, x2 and u(x3).
%lbd = [-1.0e10,-1.0e10,U(1)];
% lbd = [-1.0e10,-1.0e10,0];
lbd = [-1.0e10,-1.0e10,-1.0e10];
% upper bounds for variables x1, x2 and x3
%ubd = [1.0e10,1.0e10,U(2)];
%ubd = [1.0e10,1.0e10,U(2)-U(1)]; 
ubd = [1.0e10,1.0e10,1.0e10];

return