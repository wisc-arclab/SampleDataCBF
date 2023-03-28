% Simulation of the 2-dimension model using the sample-data CBF
clear

%% set the options of the simulation

% enable different sections of this main script
plot_flag = 1; % 1 for plotting the trajectory; 0 for no plot
clf_flag = 1; % 1 for enabling the simulation with nominal controller/ CLF controller; 0 for not
cbf_flag = 1; % 1 for enabling continuos CBF simulation; 0 for not

% set the format of the reach-tube, can be 'poly' for polytope or 'int' for interval
rtube_format = 'int';

% set the uncertainty
measure_uncertain = 1; % 1 for adding the uncertain measurement, 0 for not
if measure_uncertain
    state_eps = 0.1; % bound of the uncertain measurement
    % ''state_uncertainty_type'' options: 'max' - add positive maximum uncertainty to state (state_eps);
    %                                     'min' - add negative maximum uncertainty to state (-state_eps);
    %                                     'rand' - add bounded random uncertainty to state (between -state_eps and state_eps);
    state_uncertainty_type = 'max';
    if strcmp(state_uncertainty_type, 'rand')
        seed = 1; % seed for the 'rand' function
    end
end
actuator_uncertain = 1; % 1 for adding the uncertain actuator, 0 for not
if actuator_uncertain
    input_eps = 0.1; % bound of the uncertain actuator
    % ''state_uncertainty_type'' options: 'max' - add positive maximum uncertainty to input (input_eps);
    %                                     'min' - add negative maximum uncertainty to input (-input_eps);
    %                                     'rand' - add bounded random uncertainty to input (between -input_eps and input_eps);
    input_uncertainty_type = 'rand';
    if strcmp(input_uncertainty_type, 'rand')
        seed = 1; % seed for the 'rand' function
    end
end

% set the simulation time and sampling time
t_end = 2; % total simulation time length
dt = 0.01; % sampling time for the digital system
sim_dt = 0.01; % sampling time for solving the dynamic ODEs

% set the name of the example
param.example = 'spring'; % Current options are 'jankovic' or 'spring'

%% Load example variables and functions

if strcmp(param.example,'jankovic')
    Jankovic_example;
elseif strcmp(param.example,'spring')
    mass_spring_damper_example;
else
    error([param.example,' is not a valid example']);
end

% if consider uncertain actuator, shrink the input set U
if actuator_uncertain == 1
    U_qp = [U(:,1)+input_eps, U(:,2)-input_eps];
else
    U_qp = U;
end

%% setup quadprog & sparsePOP options

qp_options = mskoptimset(''); % get default options
qp_options = mskoptimset(qp_options,'Display','off');

sparsePOPopt.relaxOrder = 3; % relaxOrder has to be >= 3 for the 'jankovic' example and >=2 for the 'spring' example
sparsePOPopt.SDPsolver = 'mosek';
sparsePOPopt.printLevel = [0,0];
sparsePOPopt.SDPsolverOutFile = -1;
sparsePOPopt.mex = 1;
if sparsePOPopt.mex == 0
    sparsePOPopt.convert2 = 0; % This param is added by YZ to turn off convert2 (need a specific sparsePOP format). Set 0 to skip convert2, 1 not skip
end
sparsePOPopt = defaultParameter(sparsePOPopt);

%% setup the system for CORA reach

Reachparams.tFinal = dt;
Reachparams.U = zonotope(interval(U(1),U(2)));

Reachoptions.alg = 'lin';
Reachoptions.timeStep = dt;
Reachoptions.tensorOrder = 2;
Reachoptions.taylorTerms = 6;
Reachoptions.zonotopeOrder = 6;

Reachparams.R0 = zonotope(interval(x_ini,x_ini));
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
u_sim_sparse = zeros(ceil(t_end/dt),1);
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
    %tic;
    [phi,info] = sparseRange(x_k,U,param,sys,Reachoptions,sparsePOPopt,rtube_format);
    %toc;
    %phi = 0;

    % Solve the CLF-CBF-QP to get current input
    if strcmp(param.example,'jankovic')
        f = [0; 0; 0];
        A = [dVx(x_k)*gx(x_k), -1, 0; -dhx(x_k)*gx(x_k), 0, -1];
        b = [-dVx(x_k)*fx(x_k)-alpha(x_k); dhx(x_k)*fx(x_k)+param.gamma*hx(x_k)+phi];
        [res,~,exitflag,~] = quadprog(H_CLF_CBF,f,A,b,[],[],[U_qp(1);-inf;-inf],[U_qp(2);inf;inf],[],qp_options);
    elseif strcmp(param.example,'spring')
        u_ref = -1.5*x_k(1)-1.5*x_k(2);
        f = [-u_ref;0];
        A = [-dhx(x_k)*gx(x_k),-1];
        b = [dhx(x_k)*fx(x_k)+param.gamma*hx(x_k)+phi];
        [res,~,exitflag,~] = quadprog(H_CLF_CBF,f,A,b,[],[],[U_qp(1);-inf],[U_qp(2);inf],[],qp_options);
    else
        error([param.example,' is not a valid example']);
    end

    u_k = res(1);
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
    
    % uncomment below to see when "The problem is likely to be either
    % primal or dual infeasible" (see mosek documentation for quadprog)
    %{
    if exitflag < 0
        disp(['quadprog likely not feasible at time ',num2str(k)]);
        u_k
    end
    %}
    
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
    u_sim_sparse(i) = u_k;
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
    u_sim_cbf = zeros(ceil(t_end/dt),1);
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
        if strcmp(param.example,'jankovic')
            f = [0; 0; 0];
            A = [dVx(x_k)*gx(x_k), -1, 0; -dhx(x_k)*gx(x_k), 0, -1];
            b = [-dVx(x_k)*fx(x_k)-alpha(x_k); dhx(x_k)*fx(x_k)+param.gamma*hx(x_k)];
            [res,~,exitflag,~] = quadprog(H_CLF_CBF,f,A,b,[],[],[U_qp(1);-inf;-inf],[U_qp(2);inf;inf],[],qp_options);
        elseif strcmp(param.example,'spring')
            u_ref = -1.5*x_k(1)-1.5*x_k(2);
            f = [-u_ref;0];
            A = [-dhx(x_k)*gx(x_k),-1];
            b = [dhx(x_k)*fx(x_k)+param.gamma*hx(x_k)];
            [res,~,exitflag,~] = quadprog(H_CLF_CBF,f,A,b,[],[],[U_qp(1);-inf],[U_qp(2);inf],[],qp_options);
        else
            error([param.example,' is not a valid example']);
        end
        
        u_k = res(1);
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
        
        % uncomment below to see when "The problem is likely to be either
        % primal or dual infeasible" (see mosek documentation for quadprog)
        %{
    if exitflag < 0
        disp(['quadprog likely not feasible at time ',num2str(k)]);
        u_k
    end
        %}
        
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
        u_sim_cbf(i) = u_k;
        i = i+1;
    end
    
end
%% simulation with only CLF constraint or nominal controller

if clf_flag > 0
    x_k = x_ini; % store the 'estimated' state
    x_r = x_ini; % store the true state
    t_sim_clf = zeros(1+ceil(t_end/sim_dt),1);
    x_sim_clf = ones(1+ceil(t_end/sim_dt),1)*x_r';
    u_sim_clf = zeros(ceil(t_end/dt),1);
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
        
        % Solve the CLF-QP to get current input
        if strcmp(param.example,'jankovic')
            f = [0;0];
            A = [dVx(x_k)*gx(x_k), -1];
            b = [-dVx(x_k)*fx(x_k)-alpha(x_k)];
            res = quadprog(H_CLF,f,A,b,[],[],[U_qp(1);-inf],[U_qp(2);inf],[],qp_options);
            u_k = res(1);
        elseif strcmp(param.example,'spring')
            u_k = -1.5*x_k(1)-1.5*x_k(2);
            if u_k > U_qp(2)
                u_k = U_qp(2);
            elseif u_k < U_qp(1)
                u_k = U_qp(1);
            end
            tspan = k:sim_dt:k+dt;
        else
            error([param.example,' is not a valid example']);
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
        u_sim_clf(i) = u_k;
        i = i+1;
    end
end
%% Plot the constraint and trajectory
if plot_flag > 0
    if strcmp(param.example,'jankovic')
        plot(x_sim_sparse(:,1),x_sim_sparse(:,2),'Color','#0072BD','LineWidth',2);
        hold on
        if clf_flag > 0
            plot(x_sim_clf(:,1),x_sim_clf(:,2),'Color','#77AC30','LineStyle','--','LineWidth',2);
        end
        if cbf_flag > 0
            plot(x_sim_cbf(:,1),x_sim_cbf(:,2),'Color','#7E2F8E','LineStyle','-.','LineWidth',2);
        end
        x2 = -2:0.01:2;
        x1 = 1+param.q*x2.^2;
        plot(x1,x2,'r','LineWidth',2);
        xlabel('x_1');
        ylabel('x_2');
        set(gca,'FontSize', 12);
        legend('SDCBF-CLF','CBF-CLF','CLF only');
    elseif strcmp(param.example,'spring')
        %%%% Plot time and x1
        t1 = linspace(0,t_end,length(x_sim_sparse(:,1)));
        x1_top = param.x1_shift + sqrt(2*param.P_max/param.k);
        x1_bot = param.x1_shift - sqrt(2*param.P_max/param.k);
        plot(t1,x_sim_sparse(:,1),'Color','#0072BD','LineWidth',2);
        hold on;
        if clf_flag > 0
            t2 = linspace(0,t_end,length(x_sim_clf(:,1)));
            plot(t2,x_sim_clf(:,1),'Color','#77AC30','LineStyle','--','LineWidth',2);
        end
        if cbf_flag > 0
            t3 = linspace(0,t_end,length(x_sim_cbf(:,1)));
            plot(t3,x_sim_cbf(:,1),'Color','#7E2F8E','LineStyle','-.','LineWidth',2);
        end
        plot(t1,x1_top*ones(size(t1)),'r','LineWidth',2);
        plot(t1,x1_bot*ones(size(t1)),'r','LineWidth',2);
        
        xlabel('time');
        ylabel('x_1');
        set(gca,'FontSize', 12);
        legend('SDCBF + nominal','CBF + nominal','nominal controller');
        %legend('traj with initreach and CORA range bounding','traj with CLF controller','initial constraint');
%         gtext('Initial constraint','Color','red','FontSize',12);
%         gtext('CBF','Color','blue','FontSize',12);
%         gtext('Nominal control','Color','green','FontSize',12);
        
        %%%% Plot x1 and x2
        figure
        syms x1 x2
        B(x1,x2) = -param.k*x2*(x1-param.x1_shift) + param.gamma0*(-(param.k/2)*(x1 - param.x1_shift)^2 + param.P_max);
        
        hold on
        plot(x_sim_sparse(:,1),x_sim_sparse(:,2),'Color','#0072BD','LineWidth',2);
        if clf_flag > 0
            plot(x_sim_clf(:,1),x_sim_clf(:,2),'Color','#77AC30','LineStyle','--','LineWidth',2);
        end
        if cbf_flag > 0
            plot(x_sim_cbf(:,1),x_sim_cbf(:,2),'Color','#7E2F8E','LineStyle','-.','LineWidth',2);
        end
        Bc = fcontour(B,'LineColor','#EDB120');
        Bc.LevelList = [0 0];
        Bc.LineWidth = 2;
        
        plot([x1_top,x1_top],[-5,5],'r','LineWidth',2);
        plot([x1_bot,x1_bot],[-5,5],'r','LineWidth',2);
        axis([-0.5 1.5 -3 1.5]);

        xlabel('x_1');
        ylabel('x_2');
        set(gca,'FontSize', 12);
        legend('SDCBF + nominal','CBF + nominal','nominal controller');
        %legend('traj with initreach and CORA range bounding','traj with CLF controller','new constraint','original contraint');
%         gtext('Initial constraint','Color','red','FontSize',12);
%         gtext('CBF','Color','blue','FontSize',12);
%         gtext('Nominal control','Color','green','FontSize',12);
%         gtext('New constraint','Color','#09C2A9','FontSize',12);
    else
        error([param.example,' is not a valid example']);
    end
end
