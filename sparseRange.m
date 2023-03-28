function [phi,info] = sparseRange(x_k,U,param,sys,options,sparse_opt,rtube_format)
%% compute the reach-tube using CORA
t1=tic;
% compute reach tube for time interval [0, dt]
Rfirst = initReach(sys,options.R0,options);
info.reach_time = toc(t1);
%% form the polynomial range bounding problem and solve by sparsePOP
switch rtube_format
    case 'poly'
        % convert zonotopes to polytopes
        Poly = polytope(Rfirst.ti{1,1});       
        
        % Polytope set {x: Ax<=b}.
        A = Poly.P.A;
        b = Poly.P.b;
        P1 = Polyhedron(A, b);
        P2 = P1.minHRep();
        A = P2.A;
        b = P2.b;
        
        % form the sparsePOP problem
        [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_poly(A,b,x_k,U,param);
    case 'int'
        D = interval(Rfirst.ti{1,1});
        % form the sparsePOP problem
        [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_int(x_k,U,param,D);
    otherwise
        error('rtube_format shoule be either "poly" or "int".');
end

% solve the sparsePOP problem
t2=tic;
[~,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,sparse_opt);
info.POPsolving_time = toc(t2);
phi = SDPobjValue;
info.sparsePOP.elapsedTime = elapsedTime;
info.sparsePOP.SDPsolverInfo = SDPsolverInfo;
info.sparsePOP.SDPinfo = SDPinfo;

end