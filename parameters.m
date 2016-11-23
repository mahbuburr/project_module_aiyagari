%% set parameters
% consumers
beta = 0.99; % discount factor
sigma = 2; % risk aversion
mu = 0.3; % replacement rate of unemployed
k_min = 0.5; % borrowing constraint as share of capital stock
ind_no = 5000; % number of individuals simulated
T = 5000; % number of periods simulated

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
delta = 0.025; % depreciation rate
alpha = 0.33; % output elasticity of capital
delta_a=0.01;    % (1-delta_a) is the productivity level in a bad state, 
                 % and (1+delta_a) is the productivity level in a good state
z_g = 1+delta_a; % productivity in good aggreagate state
z_b = 1-delta_a; % productivity in bad aggregate state
z = [z_b; z_g];

% transition probabilities
% L_target = 0.9;
% PI_UE = 0.3; % chance of getting employed
% PI_EU = PI_UE*(1-L_target)/L_target; % chance of getting unemployed

U_b=0.1;        % unemployment rate in a bad aggregate state
L_b=(1-U_b);   % employment rate in a bad aggregate state
U_g=0.04;       % unemployment rate in a good aggregate state
L_g=(1-U_g);   % employment rate in a good aggregate state
L = [L_b; L_g]; % vector of states

%% setup
% build transition matrix of agents
% PI = [1-PI_UE,PI_UE;PI_EU,1-PI_EU]; % [UU,UE;EU;EE]
% L = PI(1,2)/(PI(1,2)+PI(2,1)); % calculate steady state supply of labour
% tau = mu*(1-L)/L; % tax rate

% define some useful functions
muc = @(c) c.^(-sigma); % marginal utility of consumption
muc_inv = @(muc) muc.^(-1/sigma); % inverse of marginal utility of consumption
w = @(K,L,z) z.*(1-alpha).*K.^alpha.*L.^(-alpha); % wage
r = @(K,L,z) z.*alpha.*K.^(alpha-1).*L.^(1-alpha); % rental rate of capital
K = @(L,r,z) L*(z*alpha/r).^(1/(1-alpha)); % capital stock necessary to yield return r
Y = @(K,L,z) z*K.^alpha*L^(1-alpha); % output
C = @(K) Y(K,L,z)-delta*K; % consumption
tau = @(L) mu*(1-L)/L; % tax rate

% define grid for individual capital on which to solve
K_rep_g = K(L_g,1/beta-1+delta, z_g);% capital stock of representative agent useful for comparison good aggregate state
K_rep_b = K(L_b,1/beta-1+delta, z_b);

grid_k_no = 100; % number of grid points for agents' capital holdings
grid_k = linspace(k_min*K_rep_b,3*K_rep_g,grid_k_no); % grid for agents' policy function

grid_dist_no = 1000;
grid_dist = linspace(grid_k(1),grid_k(end),grid_dist_no); % grid for distribution of agents - use finer grid

% Grid for mean of capital

km_min=0.75*K_rep_b;                           % minimum grid-value of the mean of 
                                     % capital distribution, km 
km_max=1.25*K_rep_g;                           % maximum grid value of km
ngridkm=4;                           % number of grid points for km 
km=linspace(km_min,km_max,ngridkm)'; % generate a grid of ngridkm points on 
                                     % [km_min,km_max] interval 

% useful matrices
mat_k = repmat(grid_k',1,2,2); % replicate grid for unemployed and employed
mat_income = @(K,L,z) w(K,L,z).*repmat([mu,1-tau(L)],grid_k_no,1); % matrix with income of each agent