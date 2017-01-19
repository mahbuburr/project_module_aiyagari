
%% set parameters
% consumers
par.beta = 0.99; % discount factor
par.sigma = 1; % risk aversion
par.mu = 0.15; % replacement rate of unemployed
par.k_min = 0; % borrowing constraint as share of capital stock

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
par.delta = 0.025; % depreciation rate
par.alpha = 0.36; % output elasticity of capital
par.z = 1; % productivity

% transition probabilities
par.L_target = 0.9;
%par.PI_EU = 0.4*(1-0.9)/0.9; % setting the layoff probability to the baseline case
par.PI_UE = 0.4; % chance of getting employed




%% solve for general equilibrium
method.HH = 'FPend'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
method.agg = 'bisection'; % 'bisection' to use bisection method, gradual updating otherwise