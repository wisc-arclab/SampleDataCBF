function [objPoly,ineqPolySys,lbd,ubd] = example112(nDim);

%
% Linearly dependent coefficient matrix, infeasible
%
% <Input> 
% nDim: The dimension of the function
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',3201);
objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 1;
objPoly.noTerms  = nDim;
objPoly.supports = speye(nDim,nDim);
objPoly.coef     = ones(nDim,1);

for i=1:nDim
ineqPolySys{i}.typeCone = -1;
ineqPolySys{i}.sizeCone = 1;
ineqPolySys{i}.dimVar   = nDim;
ineqPolySys{i}.degree   = 1;
ineqPolySys{i}.noTerms  = nDim+1;
ineqPolySys{i}.supports = [sparse(1,nDim);speye(nDim,nDim)];
ineqPolySys{i}.coef     = ones(nDim+1,1);
end
ineqPolySys{nDim}.coef(1,1) = 10;

lbd = -1.0e+10*ones(1,nDim);
ubd  =  1.0e+10*ones(1,nDim);
return;

