function [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_poly_DI(A,b,x_k,U,param)
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

%%% polytope bounds
% ineqPolySys (Ax<=b)                        
ineqPolySys = cell(size(b,1),1);       
for i = 1:size(b,1)
    ineqPolySys{i}.typeCone = 1;
    ineqPolySys{i}.dimVar   = 9;
    ineqPolySys{i}.degree   = 1;
%    ineqPolySys{i}.noTerms  = 7;
%     ineqPolySys{i}.supports = [0,0,0,0,0,0,0,0,0; 1,0,0,0,0,0,0,0,0; 0,1,0,0,0,0,0,0,0; 0,0,1,0,0,0,0,0,0; 0,0,0,1,0,0,0,0,0; 0,0,0,0,1,0,0,0,0; 0,0,0,0,0,1,0,0,0];
%     ineqPolySys{i}.coef     = [b(i); -A(i,1); -A(i,2); -A(i,3); -A(i,4); -A(i,5); -A(i,6)];
    ineqPolySys{i}.noTerms  = 1;
    ineqPolySys{i}.supports = [0,0,0,0,0,0,0,0,0];
    ineqPolySys{i}.coef     = [b(i)];
    if A(i,1) ~= 0
        ineqPolySys{i}.noTerms = ineqPolySys{i}.noTerms + 1;
        ineqPolySys{i}.supports = [ineqPolySys{i}.supports; 1,0,0,0,0,0,0,0,0];
        ineqPolySys{i}.coef = [ineqPolySys{i}.coef;-A(i,1)];
    end
    if A(i,2) ~= 0
        ineqPolySys{i}.noTerms = ineqPolySys{i}.noTerms + 1;
        ineqPolySys{i}.supports = [ineqPolySys{i}.supports; 0,1,0,0,0,0,0,0,0];
        ineqPolySys{i}.coef = [ineqPolySys{i}.coef;-A(i,2)];
    end
    if A(i,3) ~= 0
        ineqPolySys{i}.noTerms = ineqPolySys{i}.noTerms + 1;
        ineqPolySys{i}.supports = [ineqPolySys{i}.supports; 0,0,1,0,0,0,0,0,0];
        ineqPolySys{i}.coef = [ineqPolySys{i}.coef;-A(i,3)];
    end
    if A(i,4) ~= 0
        ineqPolySys{i}.noTerms = ineqPolySys{i}.noTerms + 1;
        ineqPolySys{i}.supports = [ineqPolySys{i}.supports; 0,0,0,1,0,0,0,0,0];
        ineqPolySys{i}.coef = [ineqPolySys{i}.coef;-A(i,4)];
    end
    if A(i,5) ~= 0
        ineqPolySys{i}.noTerms = ineqPolySys{i}.noTerms + 1;
        ineqPolySys{i}.supports = [ineqPolySys{i}.supports; 0,0,0,0,1,0,0,0,0];
        ineqPolySys{i}.coef = [ineqPolySys{i}.coef;-A(i,5)];
    end
    if A(i,6) ~= 0
        ineqPolySys{i}.noTerms = ineqPolySys{i}.noTerms + 1;
        ineqPolySys{i}.supports = [ineqPolySys{i}.supports; 0,0,0,0,0,1,0,0,0];
        ineqPolySys{i}.coef = [ineqPolySys{i}.coef;-A(i,6)];
    end
end
% U(:,1)<=u<=U(:,2)
for i = size(b,1)+1:size(b,1)+6
    ineqPolySys{i}.typeCone = 1;
    ineqPolySys{i}.dimVar   = 9;
    ineqPolySys{i}.degree   = 1;
    ineqPolySys{i}.noTerms  = 2;
    ineqPolySys{i}.supports = zeros(2,9);
    ineqPolySys{i}.supports(2,ceil((i-size(b,1))/2)+6) = 1;
    if mod(i-size(b,1),2)~=0
        ineqPolySys{i}.coef     = [-U(ceil((i-size(b,1))/2),1); 1];
    else
        ineqPolySys{i}.coef     = [U(ceil((i-size(b,1))/2),2); -1];
    end

end

% lower bounds for variables [X, u].
lbd = [-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10,-1.0e10];
% upper bounds for variables [X, u]
ubd = [1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10];  
        
return