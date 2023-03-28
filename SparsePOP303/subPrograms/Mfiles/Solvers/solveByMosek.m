function [x, y, SDPobjValue, SDPsolverInfo] = solveByMosek(A, b, c, K, param)
%%% YZ - 09082021: This function is not complete yet, the optimal value is
%%% returned as SDPobjValue. However, calculation for x,y,SDPsolverInfo needs further modifications.

    if ~isfield(K,'f')
       K.f = 0; 
    end
    if ~isfield(K,'l')
       K.l = 0; 
    end
    if ~isfield(K,'q')
       K.q = [];
    end
    if ~isfield(K,'s')
       K.s = [];
    end
    
    % Mosek Parameters
    % Set log level (integer parameter)
    param1.MSK_IPAR_LOG = 0;
    param1.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = param.SDPsolverEpsilon;
    param1.MSK_IPAR_WRITE_BAS_HEAD = 'MSK_OFF';
    
    % Convert sedumi problem to mosek problem
    prob = convert_sedumi2mosek(A,b,c,K);
    
    % Solve
    [Info, res] = mosekopt('minimize echo(0)',prob,param1);
    x = []; % x can not be recovered from the solution
    y = res.sol.itr.y;
    SDPsolverInfo.info = Info;
    SDPsolverInfo.res = res;
    
    if ~isfield(SDPsolverInfo, 'pinf')
        SDPsolverInfo.pinf = 0;
    end
    if ~isfield(SDPsolverInfo, 'dinf')
        SDPsolverInfo.dinf = 0;
    end
    if ~isfield(SDPsolverInfo, 'numerr')
        SDPsolverInfo.numerr = 0;
    end
    
    % The following two fomulas are equivalent to get the objective value
    SDPobjValue = b'*y;
    %SDPobjValue = res.sol.itr.pobjval;

end