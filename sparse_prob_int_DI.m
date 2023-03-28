function [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_int_DI(x_k,U,param,D)
% convert the polynomial problem with interval constraints into sparsePOP format

% objPoly 
objPoly.typeCone = 1;
objPoly.dimVar   = 9;
objPoly.degree   = 2;

% Note that the objective polynomial is changed to the equivalent one with positive variables of lower 
% bound 0 and free variables
objPoly.noTerms  = 8;
objPoly.supports = [0,0,0,0,0,0,0,0,0; 0,0,0,0,2,0,0,0,0; 0,0,0,0,1,0,0,0,0; 0,1,0,0,1,0,0,0,0; 0,2,0,0,0,0,0,0,0; 0,1,0,0,0,0,0,0,0; 0,1,0,0,0,0,0,1,0; 0,0,0,0,0,0,0,1,0];
objPoly.coef     = [-param.gamma1*param.gamma2*param.ybar*param.yunderline+2*x_k(5)^2-(param.gamma1+param.gamma2)*(param.ybar+param.yunderline)*x_k(5)+2*(param.gamma1+param.gamma2)*x_k(2)*x_k(5)-param.gamma1*param.gamma2*(param.ybar-x_k(2))*(x_k(2)-param.yunderline);...
                    -2; (param.gamma1+param.gamma2)*(param.ybar+param.yunderline); -2*(param.gamma1+param.gamma2);...
                    -(param.gamma1*param.gamma2); param.gamma1*param.gamma2*(param.ybar+param.yunderline); -2; 2*x_k(2)];

%%% for faster sparsePOP conversion, we convert the interval set into
%%% polytope constraints here

% ineqPolySys (x_inf<=x<=x_sup)                        
ineqPolySys = cell(18);       
for i=1:18
    ineqPolySys{i}.typeCone = 1;
    ineqPolySys{i}.dimVar   = 9;
    ineqPolySys{i}.degree   = 1;
    ineqPolySys{i}.noTerms  = 2;
    ineqPolySys{i}.supports = zeros(2,9);
    ineqPolySys{i}.supports(2,ceil(i/2)) = 1;
end
for i=1:6
    ineqPolySys{2*i-1}.coef     = [-D.inf(i); 1];
    ineqPolySys{2*i}.coef     = [D.sup(i); -1];
end
for i=7:9
    ineqPolySys{2*i-1}.coef     = [-U(i-6,1); 1];
    ineqPolySys{2*i}.coef     = [U(i-6,2); -1];
end

% lower bounds for variables [X, u].
lbd = [-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10];
% upper bounds for variables [X, u]
ubd = [1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10];  
        
return