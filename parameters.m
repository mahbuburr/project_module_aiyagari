%% Set parameters
% consumers
par.beta = 0.99; % discount factor
par.sigma = 2; % risk aversion
par.mu = 0.15; % replacement rate of unemployed
par.k_min = 1e-15; % borrowing constraint as share of capital stock

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
par.delta = 0.025; % depreciation rate
par.alpha = 0.36; % output elasticity of capital
par.z = 1; % productivity

% transition probabilities
%par.L_target = 0.9;
par.PI_UE = 0.764705882; % chance of getting employed
par.PI_EU = 0.08496732; % fixed layoff probability at baseline level

% solution methods
method.HH = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'simulation'; % 'simulation' for simulation, 'histogram' for histogram method
method.agg = 'bisectio'; % 'bisection' to use bisection method, gradual updating otherwise
