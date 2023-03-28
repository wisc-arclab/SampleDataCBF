function [phi,info] = sparseRange_DI(x_k,U,param,sys,options,sparse_opt,rtube_format)
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
        switch param.hxorder
            case 1
                c     = [0; -param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2);0];
                a     = A;
                blc   = [];
                buc   = b;
                blx   = [];
                bux   = []';
                t2=tic;
                res1 = msklpopt(c,a,blc,buc,blx,bux,[],'minimize echo(0)');
                res2 = msklpopt(-c,a,blc,buc,blx,bux,[],'minimize echo(0)');
                info.POPsolving_time = toc(t2);
                info.LP = {res1, res2};
                phi = [res1.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(5) + param.gamma1*param.gamma2*x_k(2);...
                       res2.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(5) - param.gamma1*param.gamma2*x_k(2)];
            case 2
                % form the sparsePOP problem
                [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_poly_DI(A,b,x_k,U,param);
            case 0
                c{1}     = [-param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2); 0; 0];
                c{2}     = [0; -param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2);0];
                c{3}     = [0; 0; -param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2)];
                a     = A;
                blc   = [];
                buc   = b;
                blx   = [];
                bux   = []';
                t2=tic;
                res{1} = msklpopt(c{1},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{2} = msklpopt(-c{1},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{3} = msklpopt(c{2},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{4} = msklpopt(-c{2},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{5} = msklpopt(c{3},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{6} = msklpopt(-c{3},a,blc,buc,blx,bux,[],'minimize echo(0)');
                info.POPsolving_time = toc(t2);
                info.LP = res;
                phi = [res{1}.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(4) + param.gamma1*param.gamma2*x_k(1);...
                       res{2}.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(4) - param.gamma1*param.gamma2*x_k(1);...
                       res{3}.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(5) + param.gamma1*param.gamma2*x_k(2);...
                       res{4}.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(5) - param.gamma1*param.gamma2*x_k(2);...
                       res{5}.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(6) + param.gamma1*param.gamma2*x_k(3);...
                       res{6}.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(6) - param.gamma1*param.gamma2*x_k(3)];
            otherwise
                error('hxorder shoule be either 0, 1 or 2.');
        end
    case 'int'
        D = interval(Rfirst.ti{1,1});
        % form the sparsePOP problem
        switch param.hxorder
            case 1
                c     = [0; -param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2);0];
                a     = eye(6);
                blc   = D.inf;
                buc   = D.sup;
                blx   = [];
                bux   = []';
                t2=tic;
                res1 = msklpopt(c,a,blc,buc,blx,bux,[],'minimize echo(0)');
                res2 = msklpopt(-c,a,blc,buc,blx,bux,[],'minimize echo(0)');
                info.POPsolving_time = toc(t2);
                info.LP = {res1, res2};
                phi = [res1.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(5) + param.gamma1*param.gamma2*x_k(2);...
                       res2.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(5) - param.gamma1*param.gamma2*x_k(2)];
            case 2
                [objPoly,ineqPolySys,lbd,ubd] = sparse_prob_int_DI(x_k,U,param,D);
            case 0
                c{1}     = [-param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2); 0; 0];
                c{2}     = [0; -param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2);0];
                c{3}     = [0; 0; -param.gamma1*param.gamma2; 0; 0; -(param.gamma1+param.gamma2)];
                a     = eye(6);
                blc   = D.inf;
                buc   = D.sup;
                blx   = [];
                bux   = []';
                t2=tic;
                res{1} = msklpopt(c{1},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{2} = msklpopt(-c{1},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{3} = msklpopt(c{2},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{4} = msklpopt(-c{2},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{5} = msklpopt(c{3},a,blc,buc,blx,bux,[],'minimize echo(0)');
                res{6} = msklpopt(-c{3},a,blc,buc,blx,bux,[],'minimize echo(0)');
                info.POPsolving_time = toc(t2);
                info.LP = res;
                phi = [res{1}.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(4) + param.gamma1*param.gamma2*x_k(1);...
                       res{2}.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(4) - param.gamma1*param.gamma2*x_k(1);...
                       res{3}.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(5) + param.gamma1*param.gamma2*x_k(2);...
                       res{4}.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(5) - param.gamma1*param.gamma2*x_k(2);...
                       res{5}.sol.itr.pobjval + (param.gamma1+param.gamma2)*x_k(6) + param.gamma1*param.gamma2*x_k(3);...
                       res{6}.sol.itr.pobjval - (param.gamma1+param.gamma2)*x_k(6) - param.gamma1*param.gamma2*x_k(3)];
            otherwise
                error('hxorder shoule be either 0, 1 or 2.');
        end
    otherwise
        error('rtube_format shoule be either "poly" or "int".');
end

if param.hxorder == 2
    % solve the sparsePOP problem
    t2=tic;
    [~,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,sparse_opt);
    info.POPsolving_time = toc(t2);
    phi = SDPobjValue;
    info.sparsePOP.elapsedTime = elapsedTime;
    info.sparsePOP.SDPsolverInfo = SDPsolverInfo;
    info.sparsePOP.SDPinfo = SDPinfo;
end

end