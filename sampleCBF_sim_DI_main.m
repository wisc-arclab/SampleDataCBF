% Simulation of the quadcopter model using the sample-data CBF

% The linearized double integrator quadcopter model
% 6-DIM Model
% state - X = [x; y; z; xdot; ydot; zdot]
% control input - u = [xddot; yddot; zddot]
% dynamics - X' = [xdot; ydot; zdot; xddot; yddot; zddot]

clc;clear all;

%% set the options of the simulation

% enable different sections of this main script
plot_flag = 1; % 1 for plotting the trajectory; 0 for no plot
clf_flag = 1; % 1 for enabling the simulation with nominal controller; 0 for not
cbf_flag = 1; % 1 for enabling continuos CBF simulation; 0 for not

% set the format of the reach-tube, can be 'poly' for polytope or 'int' for interval
rtube_format = 'int';

% set the uncertainty
measure_uncertain = 1; % 1 for adding the uncertain measurement, 0 for not
if measure_uncertain
    state_eps = 0.02; % bound of the uncertain measurement
    % ''state_uncertainty_type'' options: 'max' - add positive maximum uncertainty to state (state_eps);
    %                                     'min' - add negative maximum uncertainty to state (-state_eps);
    %                                     'rand' - add bounded random uncertainty to state (between -state_eps and state_eps);
    state_uncertainty_type = 'min';
    if strcmp(state_uncertainty_type, 'rand')
        seed = 1; % seed for the 'rand' function
    end
end
actuator_uncertain = 0; % 1 for adding the uncertain actuator, 0 for not
if actuator_uncertain
    input_eps = 0.01; % bound of the uncertain actuator
    % ''state_uncertainty_type'' options: 'max' - add positive maximum uncertainty to input (input_eps);
    %                                     'min' - add negative maximum uncertainty to input (-input_eps);
    %                                     'rand' - add bounded random uncertainty to input (between -input_eps and input_eps);
    input_uncertainty_type = 'rand';
    if strcmp(input_uncertainty_type, 'rand')
        seed = 1; % seed for the 'rand' function
    end
end

% read the nomial control sequence
traj = readmatrix('traj_fig8_20s.csv');
u_nom = traj(:,10:12);
traj_dt = traj(2,1)-traj(1,1);

% set the order of the CBF function
% param.hxorder = 0 for 6 first order constraints xunderline<=x<=xbar, yunderline<=y<=ybar, zunderline <= z <=zbar; 
% param.hxorder = 1 for first order h1(x) = ybar - x(2), h2(x) = x(2) - yunderline; 
% param.hxorder = 2 for second order h(x) = (ybar - x(2))(x(2)-yunderline)
param.hxorder = 0;

% set the simulation time and sampling time
t_end = 20; % total simulation time length (max value of t_end depends on the length of nomial control sequence from the csv file)
dt = 0.01; % sampling time for the digital system
sim_dt = 0.01; % sampling time for solving the dynamic ODEs

%% setup the linearized double integrator quadcopter model

param.g = 9.81;
param.gamma1 = 1;
param.gamma2 = 1;
param.ybar = 0.5;
param.yunderline = -0.5;
param.xbar = 0.6;
param.xunderline = -0.6;
param.zbar = 0.7;
param.zunderline = 0.0;
% Xdot = A*X + B*u = f(X)+g(X)*u.
A = [zeros(3,3), eye(3); zeros(3,6)];
B = [zeros(3,3); eye(3)];

x_ini = [0;0;0;0;0;0]; % initial state
U = [-param.g,param.g; -param.g,param.g; -param.g,param.g]/10; % constraints on inputs.
if actuator_uncertain == 1
    U_qp = [U(:,1)+input_eps, U(:,2)-input_eps];
else
    U_qp = U;
end

fx = @(x) A*[x(1),x(2),x(3),x(4),x(5),x(6)]';
gx = @(x) B;
switch param.hxorder
    case 2
        hx = @(x) (param.ybar - x(2))*(x(2)-param.yunderline);
        dhdx = @(x) [0, (param.ybar+param.yunderline)-2*x(2), 0, 0, 0, 0]; % dh/dx
        Lf_hx = @(x) (param.ybar+param.yunderline)*x(5)-2*x(2)*x(5);
        LLf_hx = @(x) -2*x(5)^2;
        Lg_hx = @(x) [0, 0, 0];
        LgLf_hx = @(x) [0,(param.ybar+param.yunderline)-2*x(2),0];
    case 1
        hx1 = @(x) param.ybar - x(2);
        dhdx1 = @(x) [0, -1, 0, 0, 0, 0]; % dh/dx
        Lf_hx1 = @(x) -x(5);
        LLf_hx1 = @(x) 0;
        Lg_hx1 = @(x) [0, 0, 0];
        LgLf_hx1 = @(x) [0,-1, 0];
        
        hx2 = @(x) x(2) - param.yunderline;
        dhdx2 = @(x) [0, 1, 0, 0, 0, 0]; % dh/dx
        Lf_hx2 = @(x) x(5);
        LLf_hx2 = @(x) 0;
        Lg_hx2 = @(x) [0, 0, 0];
        LgLf_hx2 = @(x) [0,1, 0];
    case 0
        nhx = 6;
        hx{1} = @(x) param.xbar - x(1);
        hx{2} = @(x) x(1) - param.xunderline;
        hx{3} = @(x) param.ybar - x(2);
        hx{4} = @(x) x(2) - param.yunderline;
        hx{5} = @(x) param.zbar - x(3);
        hx{6} = @(x) x(3) - param.zunderline;
        
        Lf_hx{1} = @(x) -x(4);
        Lf_hx{2} = @(x) x(4);
        Lf_hx{3} = @(x) -x(5);
        Lf_hx{4} = @(x) x(5);
        Lf_hx{5} = @(x) -x(6);
        Lf_hx{6} = @(x) x(6);
        
        for i=1:nhx
            LLf_hx{i} = @(x) 0;
        end
        
        LgLf_hx{1} = @(x) [-1, 0, 0];
        LgLf_hx{2} = @(x) [1, 0, 0];
        LgLf_hx{3} = @(x) [0, -1, 0];
        LgLf_hx{4} = @(x) [0, 1, 0];
        LgLf_hx{5} = @(x) [0, 0, -1];
        LgLf_hx{6} = @(x) [0, 0, 1];
        
    otherwise
        error('hxorder shoule be either 0, 1 or 2.');
end
% Vx = @(x) x(1)^4/4+x(2)^2/2;
% dVx = @(x) [x(1)^3, x(2)];
odefun = @(t, x, u, param) [x(4);x(5);x(6);u(1);u(2);u(3)];

%% setup quadprog & sparsePOP options

qp_options = mskoptimset(''); % get default options
qp_options = mskoptimset(qp_options,'Display','off');

sparsePOPopt.relaxOrder = 2;
sparsePOPopt.SDPsolver = 'mosek';
sparsePOPopt.printLevel = [0,0];
sparsePOPopt.SDPsolverOutFile = -1;
sparsePOPopt.mex = 1;
if sparsePOPopt.mex == 0
    sparsePOPopt.convert2 = 0; % This param is added by YZ to turn off convert2 (need a specific sparsePOP format). Set 0 to skip convert2, 1 not skip
end
sparsePOPopt = defaultParameter(sparsePOPopt);

%% setup the system for CORA reach

% dynamical model f(x)
f = @(x,u) [x(4);x(5);x(6);u(1);u(2);u(3)];
sys = nonlinearSys(f,6,3); % state dim = 6, input dim = 3.
Reachparams.tFinal = dt;
Reachparams.U = zonotope(interval(U(:,1),U(:,2)));

Reachoptions.alg = 'lin';
Reachoptions.timeStep = dt;
Reachoptions.tensorOrder = 2;
Reachoptions.taylorTerms = 6;
Reachoptions.zonotopeOrder = 6;

Reachparams.R0 = zonotope(interval(x_ini, x_ini));
Reachoptions = validateOptions(sys,'reach',Reachparams,Reachoptions);
derivatives(sys,Reachoptions);

% obtain factors for reachability analysis
for i=1:(Reachoptions.taylorTerms+1)
    %compute initial state factor
    Reachoptions.factor(i) = Reachoptions.timeStep^(i)/factorial(i);
end

%% Simulation using the sampled-data CBF constraint

x_k = x_ini; % store the 'estimated' state
x_r = x_ini; % store the true state
t_sim_sparse = zeros(1+ceil(t_end/sim_dt),1);
x_sim_sparse = ones(1+ceil(t_end/sim_dt),1)*x_r';
u_sim_sparse = zeros(ceil(t_end/dt),3);
i=1;
if measure_uncertain == 1
    if strcmp(state_uncertainty_type, 'rand')
        rng(seed);
    end
elseif actuator_uncertain == 1
    if strcmp(input_uncertainty_type, 'rand')
        rng(seed);
    end
end

sim_start=tic;
for k=0:dt:(t_end-dt)
    
    if measure_uncertain == 1
        % if there is uncertain measurement then x_k is the estimated state
        switch state_uncertainty_type
            case 'max'
                x_k = x_r + state_eps; % max uncertainty
            case 'min'
                x_k = x_r - state_eps; % - max uncertainty
            case 'rand'
                x_k = x_r + 2*state_eps*rand(size(x_r))-state_eps;
            otherwise
                error('state_uncertainty_type shoule be "max", "min" or "rand".');
        end
        Reachoptions.R0 = zonotope(interval(x_k-state_eps, x_k+state_eps));
    else
        x_k = x_r;
        Reachoptions.R0 = zonotope(interval(x_k, x_k));
    end

    % get lower bound of the CBF constraint
    %t1=tic;
    [phi,info] = sparseRange_DI(x_k,U,param,sys,Reachoptions,sparsePOPopt,rtube_format);
    %toc(t1);
    
    % Solve the CLF-CBF-QP to get current input
    switch param.hxorder
        case 2
            H = diag([1,1,1,10]);
            f = [-u_nom(floor(k/traj_dt)+1,:)'; 0];
            A = [ -LgLf_hx(x_k), -1];
            b = [LLf_hx(x_k)+(param.gamma1+param.gamma2)*Lf_hx(x_k)+param.gamma1*param.gamma2*hx(x_k)+phi];
            res = quadprog(H,f,A,b,[],[],[U_qp(:,1);-inf],[U_qp(:,2);inf],[],qp_options);
        case 1
            H = diag([1,1,1,10,10]);
            f = [-u_nom(floor(k/traj_dt)+1,:)'; 0; 0];
            A = [ -LgLf_hx1(x_k), -1, 0; -LgLf_hx2(x_k), 0, -1];
            b = [LLf_hx1(x_k)+(param.gamma1+param.gamma2)*Lf_hx1(x_k)+param.gamma1*param.gamma2*hx1(x_k)+phi(1); ...
                 LLf_hx2(x_k)+(param.gamma1+param.gamma2)*Lf_hx2(x_k)+param.gamma1*param.gamma2*hx2(x_k)+phi(2)];
            res = quadprog(H,f,A,b,[],[],[U_qp(:,1);-inf;-inf],[U_qp(:,2);inf;inf],[],qp_options);
        case 0
            H = diag([1,1,1,10*ones(1,nhx)]);
            f = [-u_nom(floor(k/traj_dt)+1,:)'; zeros(nhx,1)];
            A = zeros(nhx,nhx+3);
            b = zeros(nhx,1);
            for j = 1:nhx
                A(j,1:3) = -LgLf_hx{j}(x_k);
                A(j,j+3) = -1;
                b(j,1) = LLf_hx{j}(x_k)+(param.gamma1+param.gamma2)*Lf_hx{j}(x_k)+param.gamma1*param.gamma2*hx{j}(x_k)+phi(j);
            end
            res = quadprog(H,f,A,b,[],[],[U_qp(:,1);-inf*ones(nhx,1)],[U_qp(:,2);inf*ones(nhx,1)],[],qp_options);
        otherwise
            error('hxorder shoule be either 0, 1 or 2.');
    end
    
    u_k = res(1:3);
    if actuator_uncertain
        switch input_uncertainty_type
            case 'max'
                u_k = u_k + input_eps; % max uncertainty
            case 'min'
                u_k = u_k - input_eps; % -max uncertainty
            case 'rand'
                u_k = u_k + 2*input_eps*rand(size(u_k))-input_eps; 
            otherwise
                error('input_uncertainty_type shoule be "max", "min" or "rand".');
        end
    end
    
    % apply the input to the system for dt sec
    tspan = k:sim_dt:k+dt;
    [t,x] = ode45(@(t,x) odefun(t,x,u_k,param), tspan, x_r, odeset('RelTol',1e-3));
    
    % store the trajectory
    if size(tspan,2) == 2
        t_sim_sparse(i+1) = t(end);
        x_sim_sparse(i+1,:) = x(end,:);
    else
        t_sim_sparse((i-1)*ceil(dt/sim_dt)+2:i*ceil(dt/sim_dt)+1) = t(2:end);
        x_sim_sparse((i-1)*ceil(dt/sim_dt)+2:i*ceil(dt/sim_dt)+1,:) = x(2:end,:);
    end
    
    x_r = x(end,:)';
    u_sim_sparse(i,:) = u_k';
    i = i+1;
end
dur = toc(sim_start);
disp(['Total simulation time is ',num2str(dur),' sec.']);

%% Simulation using the continuous CBF constraint

if cbf_flag >0
    x_k = x_ini; % store the 'estimated' state
    x_r = x_ini; % store the true state
    t_sim_cbf = zeros(1+ceil(t_end/sim_dt),1);
    x_sim_cbf = ones(1+ceil(t_end/sim_dt),1)*x_r';
    u_sim_cbf = zeros(ceil(t_end/dt),3);
    i=1;
    if measure_uncertain == 1
        if strcmp(state_uncertainty_type, 'rand')
            rng(seed);
        end
    elseif actuator_uncertain == 1
        if strcmp(input_uncertainty_type, 'rand')
            rng(seed);
        end
    end
    
    for k=0:dt:(t_end-dt)
        
        if measure_uncertain == 1
            % if there is uncertain measurement then x_k is the estimated state
            switch state_uncertainty_type
                case 'max'
                    x_k = x_r + state_eps; % max uncertainty
                case 'min'
                    x_k = x_r - state_eps; % - max uncertainty
                case 'rand'
                    x_k = x_r + 2*state_eps*rand(size(x_r))-state_eps;
                otherwise
                    error('state_uncertainty_type shoule be "max", "min" or "rand".');
            end
        else
            x_k = x_r;
        end
        
        % Solve the CLF-CBF-QP to get current input
        switch param.hxorder
        case 2
            H = diag([1,1,1,10]);
            f = [-u_nom(floor(k/traj_dt)+1,:)'; 0];
            A = [ -LgLf_hx(x_k), -1];
            b = [LLf_hx(x_k)+(param.gamma1+param.gamma2)*Lf_hx(x_k)+param.gamma1*param.gamma2*hx(x_k)];
            res = quadprog(H,f,A,b,[],[],[U_qp(:,1);-inf],[U_qp(:,2);inf],[],qp_options);
        case 1
            H = diag([1,1,1,10,10]);
            f = [-u_nom(floor(k/traj_dt)+1,:)'; 0; 0];
            A = [ -LgLf_hx1(x_k), -1, 0; -LgLf_hx2(x_k), 0, -1];
            b = [LLf_hx1(x_k)+(param.gamma1+param.gamma2)*Lf_hx1(x_k)+param.gamma1*param.gamma2*hx1(x_k); ...
                 LLf_hx2(x_k)+(param.gamma1+param.gamma2)*Lf_hx2(x_k)+param.gamma1*param.gamma2*hx2(x_k)];
            res = quadprog(H,f,A,b,[],[],[U_qp(:,1);-inf;-inf],[U_qp(:,2);inf;inf],[],qp_options);
        case 0
            H = diag([1,1,1,10*ones(1,nhx)]);
            f = [-u_nom(floor(k/traj_dt)+1,:)'; zeros(nhx,1)];
            A = zeros(nhx,nhx+3);
            b = zeros(nhx,1);
            for j = 1:nhx
                A(j,1:3) = -LgLf_hx{j}(x_k);
                A(j,j+3) = -1;
                b(j,1) = LLf_hx{j}(x_k)+(param.gamma1+param.gamma2)*Lf_hx{j}(x_k)+param.gamma1*param.gamma2*hx{j}(x_k);
            end
            res = quadprog(H,f,A,b,[],[],[U_qp(:,1);-inf*ones(nhx,1)],[U_qp(:,2);inf*ones(nhx,1)],[],qp_options);
        otherwise
            error('hxorder shoule be either 0, 1 or 2.');
        end
        
        u_k = res(1:3);
        if actuator_uncertain
            switch input_uncertainty_type
                case 'max'
                    u_k = u_k + input_eps; % max uncertainty
                case 'min'
                    u_k = u_k - input_eps; % -max uncertainty
                case 'rand'
                    u_k = u_k + 2*input_eps*rand(size(u_k))-input_eps;
                otherwise
                    error('input_uncertainty_type shoule be "max", "min" or "rand".');
            end
        end
        
        % apply the input to the system for dt sec
        tspan = k:sim_dt:k+dt;
        [t,x] = ode45(@(t,x) odefun(t,x,u_k,param), tspan, x_r, odeset('RelTol',1e-3));
        
        % store the trajectory
        if size(tspan,2) == 2
            t_sim_cbf(i+1) = t(end);
            x_sim_cbf(i+1,:) = x(end,:);
        else
            t_sim_cbf((i-1)*ceil(dt/sim_dt)+2:i*ceil(dt/sim_dt)+1) = t(2:end);
            x_sim_cbf((i-1)*ceil(dt/sim_dt)+2:i*ceil(dt/sim_dt)+1,:) = x(2:end,:);
        end
        
        x_r = x(end,:)';
        u_sim_cbf(i,:) = u_k';
        i = i+1;
    end
end

%% simulation with nominal controller

if clf_flag > 0
    x_r = x_ini; % store the true state
    t_sim_clf = zeros(1+ceil(t_end/sim_dt),1);
    x_sim_clf = ones(1+ceil(t_end/sim_dt),1)*x_r';
    u_sim_clf = zeros(ceil(t_end/dt),3);
    i=1;
    
    for k=0:dt:(t_end-dt)
        u_k = u_nom(floor(k/traj_dt)+1,:)';
        for j=1:size(u_k,1)
            if u_k(j) > U_qp(j,2)
                u_k(j) = U_qp(j,2);
            elseif u_k(j) < U_qp(j,1)
                u_k(j) = U_qp(j,1);
            end
        end
        
        if actuator_uncertain
            switch input_uncertainty_type
                case 'max'
                    u_k = u_k + input_eps; % max uncertainty
                case 'min'
                    u_k = u_k - input_eps; % -max uncertainty
                case 'rand'
                    u_k = u_k + 2*input_eps*rand(size(u_k))-input_eps;
                otherwise
                    error('input_uncertainty_type shoule be "max", "min" or "rand".');
            end
        end
        tspan = k:sim_dt:k+dt;
        
        % apply the input to the system for dt sec
        [t,x] = ode45(@(t,x) odefun(t,x,u_k,param), tspan, x_r, odeset('RelTol',1e-3));
        
        % store the trajectory
        if size(tspan,2) == 2
            t_sim_clf(i+1) = t(end);
            x_sim_clf(i+1,:) = x(end,:);
        else
            t_sim_clf((i-1)*ceil(dt/sim_dt)+2:i*ceil(dt/sim_dt)+1) = t(2:end);
            x_sim_clf((i-1)*ceil(dt/sim_dt)+2:i*ceil(dt/sim_dt)+1,:) = x(2:end,:);
        end
        
        x_r = x(end,:)';
        u_sim_clf(i,:) = u_k';
        i = i+1;
    end
end

%% Plot the constraint and trajectory
if plot_flag > 0
    %%%% Plot time and x2
    t1 = linspace(0,t_end,length(x_sim_sparse(:,2)));
    plot(t1,x_sim_sparse(:,2),'Color','#0072BD','LineWidth',2);
    hold on;
    if clf_flag > 0
        t2 = linspace(0,t_end,length(x_sim_clf(:,2)));
        %plot(t2,x_sim_clf(:,2),'Color','#77AC30','LineStyle','--','LineWidth',2);
    end
    if cbf_flag > 0
        t3 = linspace(0,t_end,length(x_sim_cbf(:,2)));
        plot(t3,x_sim_cbf(:,2),'Color','#7E2F8E','LineStyle','-.','LineWidth',2);
    end
    plot(t1,param.ybar*ones(size(t1)),'r','LineWidth',2);
    plot(t1,param.yunderline*ones(size(t1)),'r','LineWidth',2);
    xlabel('time');
    ylabel('y');
    set(gca,'FontSize', 12);
    legend('SDCBF + nominal','CBF + nominal','nominal controller');
%     t1 = linspace(0,t_end,length(x_sim_sparse(:,2)));
%     plot(t1,x_sim_sparse(:,2),'Color','#0072BD','LineWidth',2);
%     hold on;
%     t3 = linspace(0,t_end,length(x_sim_cbf(:,2)));
%     plot(t3,x_sim_cbf(:,2),'Color','#7E2F8E','LineStyle','-.','LineWidth',2);
%     
%     t1 = linspace(0,t_end,length(x_sim_sparse(:,2)));
%     t3 = linspace(0,t_end,length(x_sim_cbf(:,2)));
%     plot(t1,x_sim_sparse(:,2),'Color','#EDB120','LineWidth',2);plot(t3,x_sim_cbf(:,2),'Color','#D95319','LineStyle','-.','LineWidth',2);
%     
%     t1 = linspace(0,t_end,length(x_sim_sparse(:,2)));
%     t3 = linspace(0,t_end,length(x_sim_cbf(:,2)));
%     plot(t1,x_sim_sparse(:,2),'Color','#77AC30','LineWidth',2);plot(t3,x_sim_cbf(:,2),'Color','#4DBEEE','LineStyle','-.','LineWidth',2);
    
    figure
    plot3(x_sim_sparse(:,1),x_sim_sparse(:,2),x_sim_sparse(:,3),'Color','#0072BD','Linewidth',2);
    hold on
    if clf_flag > 0
        plot3(x_sim_clf(:,1),x_sim_clf(:,2),x_sim_clf(:,3),'Color','#77AC30','Linewidth',2);
    end
    if cbf_flag > 0
        plot3(x_sim_cbf(:,1),x_sim_cbf(:,2),x_sim_cbf(:,3),'Color','#7E2F8E','Linewidth',2);
    end
    if param.hxorder == 0
        hr = [param.xunderline, param.yunderline, param.zunderline;...
            param.xbar, param.yunderline, param.zunderline;...
            param.xbar, param.ybar, param.zunderline;...
            param.xunderline, param.ybar, param.zunderline;...
            param.xunderline, param.yunderline, param.zbar;...
            param.xbar, param.yunderline, param.zbar;...
            param.xbar, param.ybar, param.zbar;...
            param.xunderline, param.ybar, param.zbar];
        
        face = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
        patch('Vertices',hr,'Faces',face,'FaceColor','r','FaceAlpha',0.3);
    else
        [xmesh,zmesh] = meshgrid(-1:2:1,-0.2:1:0.8);
        ybar = param.ybar*ones(size(xmesh));
        surf(xmesh,ybar,zmesh,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);
        yunder = param.yunderline*ones(size(xmesh));
        surf(xmesh,yunder,zmesh,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
    set(gca,'FontSize', 12);
    legend('SDCBF + nominal','CBF + nominal','nominal controller');
end