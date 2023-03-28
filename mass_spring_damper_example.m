% Mass-spring-damper example info

% Variables
x_ini = [0.5;-1]; % initial state
U = [-10,10]; % constraints on inputs. if X-ini= [-4,1], U = [-30,1]

param.k = 50;
param.gamma0 = 10;
param.gamma = 10;
param.P_max = 10;
param.mass = 1.5;
param.b = 3;
param.x1_shift = x_ini(1);

% bounds on [x1,x2,u]
param.lbd = [-.2,-3,U(1)];  % optional
param.ubd = [1.2,2,U(2)];

% Symbolic variables for check in sparse_prob.m (usually commented)
%{
syms k P_max mass b x1_shift gamma0
param.k = k;
param.P_max = P_max;
param.mass = mass;
param.b = b;
param.x1_shift = x1_shift;
param.gamma0 = gamma0;
%}

A_sys= [0,1; -param.k/param.mass,-param.b/param.mass];
B_sys = [0; 1/param.mass];

% Functions
fx = @(x) A_sys*[x(1); x(2)];
gx = @(x) B_sys;
dhx = @(x) -param.k*[x(2)+param.gamma0*(x(1)-param.x1_shift), x(1)-param.x1_shift];
hx = @(x) -param.k*x(2)*(x(1)-param.x1_shift) + param.gamma0*(-(param.k/2)*(x(1) - param.x1_shift)^2 + param.P_max);
Vx = @(x) 0.5*(param.k*x(1)^2 + param.mass*x(2)^2);
dVx = @(x) [param.k*x(1), param.mass*x(2)];
alpha = @(x) 0.5*param.mass*x(2)^2; 

% Dynamical model f(x), needed for CORA
f = @(x,u) A_sys*[x(1);x(2)] + B_sys*u;
sys = nonlinearSys(f,2,1); % state dim = 2, input dim = 1.

% H matrices for the QPs
% H_CLF_CBF = [1 0 0; 0 1 0; 0 0 1]; % Use this matrix when using clf
H_CLF_CBF = [1,0;0,1];
H_CLF = [1 0;0 8];

% ODE function
odefun = @(t, x, u, param) [0,1; -param.k/param.mass,-param.b/param.mass]*[x(1);x(2)] + ...
    [0; 1/param.mass]*u;
