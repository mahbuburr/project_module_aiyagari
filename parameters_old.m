
%% set parameters
% consumers
par.beta = 0.99; % discount factor
par.sigma = 2; % risk aversion
par.mu = 0.3; % replacement rate of unemployed
par.k_min = 0.5; % borrowing constraint as share of capital stock

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
par.delta = 0.025; % depreciation rate
par.alpha = 0.33; % output elasticity of capital
par.z = 1; % productivity

% transition probabilities
par.L_target = 0.9;
par.PI_UE = 0.3; % chance of getting employed




%% solve for general equilibrium
method.HH = 'FPend'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
method.agg = 'bisection'; % 'bisection' to use bisection method, gradual updating otherwise

%% setup
par.PI_EU = par.PI_UE*(1-par.L_target)/par.L_target; % chance of getting unemployed

% build transition matrix of agents
par.PI = [1-par.PI_UE,par.PI_UE;par.PI_EU,1-par.PI_EU]; % [UU,UE;EU;EE]
par.L = par.PI(1,2)/(par.PI(1,2)+par.PI(2,1)); % calculate steady state supply of labour
par.tau = par.mu*(1-par.L)/par.L; % tax rate

% define some useful functions
func.muc = @(c) c.^(-par.sigma); % marginal utility of consumption
func.muc_inv = @(muc) muc.^(-1/par.sigma); % inverse of marginal utility of consumption
func.w = @(K) par.z*(1-par.alpha)*K.^par.alpha*par.L^(-par.alpha); % wage
func.r = @(K) par.z*par.alpha*K.^(par.alpha-1)*par.L^(1-par.alpha); % rental rate of capital
func.K = @(r) par.L*(par.z*par.alpha/r).^(1/(1-par.alpha)); % capital stock necessary to yield return r
func.Y = @(K) par.z*K.^par.alpha*par.L^(1-par.alpha); % output
func.C = @(K) func.Y(K)-par.delta*K; % consumption
func.U = @(c) (c.^(1-par.sigma))./(1-par.sigma); % Utility
