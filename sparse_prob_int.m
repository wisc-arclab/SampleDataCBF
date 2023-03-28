function [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_int(x_k,U,param,interval)
% convert the polynomial problem with interval constraints into sparsePOP format

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

%%% for faster sparsePOP conversion, we convert the interval set into
%%% polytope constraints here

% ineqPolySys (x_inf<=x<=x_sup)                        
ineqPolySys = cell(6);       
for i=1:6
    ineqPolySys{i}.typeCone = 1;
    ineqPolySys{i}.dimVar   = 3;
    ineqPolySys{i}.degree   = 1;
    ineqPolySys{i}.noTerms  = 2;
end
ineqPolySys{1}.supports = [0,0,0; 1,0,0];
ineqPolySys{1}.coef     = [-interval.inf(1); 1];
ineqPolySys{2}.supports = [0,0,0; 1,0,0];
ineqPolySys{2}.coef     = [interval.sup(1); -1];
ineqPolySys{3}.supports = [0,0,0; 0,1,0];
ineqPolySys{3}.coef     = [-interval.inf(2); 1];
ineqPolySys{4}.supports = [0,0,0; 0,1,0];
ineqPolySys{4}.coef     = [interval.sup(2); -1];
ineqPolySys{5}.supports = [0,0,0; 0,0,1];
ineqPolySys{5}.coef     = [-U(1); 1];
ineqPolySys{6}.supports = [0,0,0; 0,0,1];
ineqPolySys{6}.coef     = [U(2); -1];

% lower bounds for variables x1, x2 and u(x3).
lbd = [-1.0e10,-1.0e10,-1.0e10];
% upper bounds for variables x1, x2 and x3
ubd = [1.0e10,1.0e10,1.0e10];  
        
return