%% set parameters
% consumers
beta = 0.99;      % discount factor
sigma = 1;        % risk aversion
mu = 0.15;        % replacement rate of unemployed
k_min = 10e-10;   % borrowing constraint as share of capital stock
ind_no = 10000;    % number of individuals simulated
T = 1100;         % number of periods simulated
e = [0,1];        % employment states
e_s = [1;2];

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
delta = 0.025;   % depreciation rate
alpha = 0.36;    % output elasticity of capital
delta_a=0.01;    % (1-delta_a) is the productivity level in a bad state, 
                 % and (1+delta_a) is the productivity level in a good state
z_g = 1+delta_a; % productivity in good aggreagate state
z_b = 1-delta_a; % productivity in bad aggregate state
z = [z_b; z_g];  % vetor of productivity states
z_s = [1;2];
ag_states_no = 2;% number of aggregate states
id_states_no = 2;% number of idiosyncratic states

% transition probabilities

U_b=0.1;        % unemployment rate in a bad aggregate state
L_b=(1-U_b);    % employment rate in a bad aggregate state
PI_UE_b = 0.4;  % chance of getting employed in bad state

U_g=0.04;       % unemployment rate in a good aggregate state
L_g=(1-U_g);    % employment rate in a good aggregate state
PI_UE_g = 2/3;  % chance of getting employed in good state

L = [L_b; L_g]; % vector of states for labour
l_bar=1/L_b;    % used for simplification

B_p = 0.5; % bad states proportion
PI_bg = 0.125; % chance of observing a bad state when at good state
%% setup

% define some useful functions
muc = @(c) c.^(-sigma);                                             % marginal utility of consumption
muc_inv = @(muc) muc.^(-1/sigma);                                   % inverse of marginal utility of consumption
w = @(K,L,z) z.*(1-alpha).*K.^alpha.*(L*l_bar).^(-alpha);           % wage
w_mat = @(K,L,z) repmat(w(K',L,z),100,1,2,2);                       % wage matrix used for aggregate uncertainty
r = @(K,L,z) z.*alpha.*K.^(alpha-1).*(L*l_bar).^(1-alpha);          % rental rate of capital
r_mat = @(K,L,z) (1-delta+repmat(r(K',L,z),100,1,2,2));             % interest rate matrix used for aggregate uncertainty
K = @(L,r,z) L*(z*alpha/r).^(1/(1-alpha));                          % capital stock necessary to yield return r
Y = @(K,L,z) z*K.^alpha*L^(1-alpha);                                % output
C = @(K) Y(K,L,z)-delta*K;                                          % consumption
tau = @(L) mu*(1-L)/(l_bar*L);                                      % tax rate

% define grid for individual capital on which to solve
grid_k_no = 100;                        % number of grid points for agents' capital holdings
% K_rep_g = K(L_g,1/beta-1+delta, z_g);   % capital stock of representative agent useful for comparison good aggregate state
% K_rep_b = K(L_b,1/beta-1+delta, z_b);   % capital stock of representative agent useful for comparison bad aggregate state
% grid_k = linspace(k_min*K_rep_b,3*K_rep_g,grid_k_no); % grid for agents' policy function
k_min=0;                      % minimum grid-value of capital
k_max=1000;                   % maximum grid-value of capital
x=linspace(0,0.5,grid_k_no)'; % generate a grid of ngridk points on [0,0.5] 
y=x.^7/max(x.^7);             % polynomial distribution of grid points 
grid_k=k_min+(k_max-k_min)*y; % transformation of grid points from [0,0.5] interval to [k_min,k_max] interval

grid_dist_no = 1000;
grid_dist = linspace(grid_k(1),grid_k(end),grid_dist_no); % grid for distribution of agents - use finer grid

% Grid for mean of capital
K_min = 30;         
K_max = 50;
% K_min=0.75*K_rep_b;  % minimum grid-value of the mean of capital distribution, K
% K_max=1.25*K_rep_g;  % maximum grid value of km
grid_K_no=4;                             % number of grid points for K 
grid_K=linspace(K_min,K_max,grid_K_no)'; % generate a grid of grid_K_no points on [K_min,K_max] interval

% useful matricies
mat_k = repmat(grid_k,1,4,2,2);                                     % replicate grid for unemployed and employed
mat_income = @(K,L,z,e) ((1-tau(L))*l_bar*e + mu*(1-e)).*w_mat(K,L,z); % matrix with income of each agent
income = @(K,L,z,e) ((1-tau(L))*l_bar*e + mu*(1-e)).*w(K',L,z);

% Forecasting
B=[0 1 0 1];    
ones4 = ones(100,4,2,2);  

%% Aggregate problem
kss=((1/beta-(1-delta))/alpha)^(1/(alpha-1));
update_B=0.3;