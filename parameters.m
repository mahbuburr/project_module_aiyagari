%% set parameters
% consumers
beta = 0.99; % discount factor
sigma = 2; % risk aversion
mu = 0.2; % replacement rate of unemployed
k_min = 0.5; % borrowing constraint as share of capital stock
ind_no = 5000; % number of individuals simulated
T = 5000; % number of periods simulated

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
delta = 0.025; % depreciation rate
alpha = 0.33; % output elasticity of capital
z = 1; % productivity

% transition probabilities
L_target = 0.9;
PI_UE = 0.3; % chance of getting employed
PI_EU = PI_UE*(1-L_target)/L_target; % chance of getting unemployed


%% setup
% build transition matrix of agents
PI = [1-PI_UE,PI_UE;PI_EU,1-PI_EU]; % [UU,UE;EU;EE]
L = PI(1,2)/(PI(1,2)+PI(2,1)); % calculate steady state supply of labour
tau = mu*(1-L)/L; % tax rate

% define some useful functions
muc = @(c) c.^(-sigma); % marginal utility of consumption
muc_inv = @(muc) muc.^(-1/sigma); % inverse of marginal utility of consumption
w = @(K) z*(1-alpha)*K.^alpha*L^(-alpha); % wage
r = @(K) z*alpha*K.^(alpha-1)*L^(1-alpha); % rental rate of capital
K = @(r) L*(z*alpha/r).^(1/(1-alpha)); % capital stock necessary to yield return r
Y = @(K) z*K.^alpha*L^(1-alpha); % output
C = @(K) Y(K)-delta*K; % consumption

% define grid for individual capital on which to solve
K_rep = K(1/beta-1+delta);% capital stock of representative agent useful for comparison

grid_k_no = 100; % number of grid points for agents' capital holdings
grid_k = linspace(k_min*K_rep,3*K_rep,grid_k_no); % grid for agents' policy function

grid_dist_no = 1000;
grid_dist = linspace(grid_k(1),grid_k(end),grid_dist_no); % grid for distribution of agents - use finer grid

% useful matrices
mat_k = [grid_k',grid_k']; % replicate grid for unemployed and employed
mat_income = @(K) w(K)*repmat([mu,1-tau],grid_k_no,1); % matrix with income of each agent