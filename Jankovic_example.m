% Jankovic example info

% This file contains variables and funtions from the example in
% "Robust control barrier functions for constrained stabilization of
% nonlinear systems.

% Variables for example
param.d = 0.6;
param.q = -1;
param.gamma = 3;
x_ini = [-2;1]; % initial state
U = [-1,1]; % constraints on inputs. if X-ini= [-4,1], U = [-30,1]

% Symbolic variables for check in sparse_prob.m (usually commented)
%{
syms d q
param.d = d;
param.q = q;
%}

% Functions
fx = @(x) [-param.d*x(1)-x(2); x(1)^3];
gx = @(x) [0 ; x(2)];
dhx = @(x) [-1, param.q*2*x(2)];
hx = @(x) param.q*x(2)^2-x(1)+1;
Vx = @(x) x(1)^4/4+x(2)^2/2;
dVx = @(x) [x(1)^3, x(2)];
alpha = @(x_k) param.d/2*(norm(x_k))^4;

% Dynamical model f(x), needed for CORA
f = @(x,u) [-param.d*x(1)-x(2); 
            x(1)^3+x(2)*u];
sys = nonlinearSys(f,2,1); % state dim = 2, input dim = 1.

% H matrices for the QPs
H_CLF_CBF = [1 0 0; 0 80 0; 0 0 8];
H_CLF = [1 0;0 8];

% ODE function
odefun = @(t, x, u, param) [-param.d*x(1)-x(2); x(1)^3+x(2)*u];
