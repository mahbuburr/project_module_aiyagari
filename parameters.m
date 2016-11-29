%% set parameters
% consumers
beta = 0.99; % discount factor
sigma = 1; % risk aversion
mu = 0.15; % replacement rate of unemployed
k_min = 10e-10; % borrowing constraint as share of capital stock
ind_no = 5000; % number of individuals simulated
T = 3000; % number of periods simulated

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
delta = 0.025; % depreciation rate
alpha = 0.36; % output elasticity of capital
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
l_bar=1/0.9;

%% setup
% build transition matrix of agents
% PI = [1-PI_UE,PI_UE;PI_EU,1-PI_EU]; % [UU,UE;EU;EE]
% L = PI(1,2)/(PI(1,2)+PI(2,1)); % calculate steady state supply of labour
% tau = mu*(1-L)/L; % tax rate

% define some useful functions
muc = @(c) c.^(-sigma); % marginal utility of consumption
muc_inv = @(muc) muc.^(-1/sigma); % inverse of marginal utility of consumption
w = @(K,L,z) z.*(1-alpha).*K.^alpha.*(L*l_bar).^(-alpha); % wage
r = @(K,L,z) z.*alpha.*K.^(alpha-1).*(L*l_bar).^(1-alpha); % rental rate of capital
K = @(L,r,z) L*(z*alpha/r).^(1/(1-alpha)); % capital stock necessary to yield return r
Y = @(K,L,z) z*K.^alpha*L^(1-alpha); % output
C = @(K) Y(K,L,z)-delta*K; % consumption
tau = @(L) mu*(1-L)/L; % tax rate

r_mat = @(K,L,z) repmat(r(K',L,z),100,1,2,2);
w_mat = @(K,L,z) repmat(w(K',L,z),100,1,2,2);

% define grid for individual capital on which to solve
K_rep_g = K(L_g,1/beta-1+delta, z_g);% capital stock of representative agent useful for comparison good aggregate state
K_rep_b = K(L_b,1/beta-1+delta, z_b);

grid_k_no = 100; % number of grid points for agents' capital holdings
% grid_k = linspace(k_min*K_rep_b,3*K_rep_g,grid_k_no); % grid for agents' policy function

k_min=0;                   % minimum grid-value of capital
k_max=1000;                % maximum grid-value of capital                % number of grid points
x=linspace(0,0.5,grid_k_no)'; % generate a grid of ngridk points on [0,0.5] 
                           % interval  
y=x.^7/max(x.^7);          % polynomial distribution of grid points, formula 
                           % (7) in the paper
grid_k=k_min+(k_max-k_min)*y;   % transformation of grid points from [0,0.5] 
                           % interval to [k_min,k_max] interval

grid_dist_no = 1000;
grid_dist = linspace(grid_k(1),grid_k(end),grid_dist_no); % grid for distribution of agents - use finer grid

% Grid for mean of capital

% km_min=0.75*K_rep_b;                           % minimum grid-value of the mean of 
km_min = 30;        % capital distribution, km 
km_max = 50;
% km_max=1.25*K_rep_g;                           % maximum grid value of km
ngridkm=4;                           % number of grid points for km 
km=linspace(km_min,km_max,ngridkm)'; % generate a grid of ngridkm points on 
                                     % [km_min,km_max] interval 

% useful matrices
mat_k = repmat(grid_k,1,4,2,2); % replicate grid for unemployed and employed
mat_income = @(K,L,z) w(K,L,z).*repmat([mu,1-tau(L)],grid_k_no,1); % matrix with income of each agent

% Forecasting
B=[0 1 0 1];
kmaux=zeros(grid_k_no,ngridkm,2,2); % for the mean of capital 
                                                   % distribution (km)
kmaux(:,:,1,1)=ones(grid_k_no,1)*km';
kmaux(:,:,1,2)=ones(grid_k_no,1)*km';
kmaux(:,:,2,1)=ones(grid_k_no,1)*km';
kmaux(:,:,2,2)=ones(grid_k_no,1)*km';

kmprime=zeros(grid_k_no,ngridkm,2,2);

kmprime(:,:,1,1)=exp(B(1)*ones(grid_k_no,ngridkm)+B(2)*log(kmaux(:,:,1,1)));
kmprime(:,:,1,2)=exp(B(1)*ones(grid_k_no,ngridkm)+B(2)*log(kmaux(:,:,1,2)));
kmprime(:,:,2,1)=exp(B(3)*ones(grid_k_no,ngridkm)+B(4)*log(kmaux(:,:,2,1)));
kmprime(:,:,2,2)=exp(B(3)*ones(grid_k_no,ngridkm)+B(4)*log(kmaux(:,:,2,2)));
kmprime=(kmprime>=km_min).*(kmprime<=km_max).*kmprime+(kmprime<km_min)*km_min+(kmprime>km_max)*km_max; % restricting km' to be in [km_min,km_max] range

nstates_ag = 2;
nstates_id = 2;
prob_bu=zeros(grid_k_no,4,nstates_ag,nstates_id); % for a bad agg. state 
                           % and unemployed idios. state in the next period
prob_be=zeros(grid_k_no,ngridkm,nstates_ag,nstates_id); % for a bad agg. state 
                           % and employed idios. state in the next period
prob_gu=zeros(grid_k_no,ngridkm,nstates_ag,nstates_id); % for a good agg. state 
                           % and unemployed idios. state in the next period
prob_ge=zeros(grid_k_no,ngridkm,nstates_ag,nstates_id); % for a good agg. state 
                           % and employed idios. state in the next period

%% Parameters for Krussel - Smith

% Matrix of transition probabilities in Den Haan, Judd, Juillard (2008)

prob=[0.525 0.35 0.03125 0.09375  
   0.038889 0.836111 0.002083 0.122917
   0.09375 0.03125 0.291667 0.583333
   0.009115 0.115885 0.024306 0.850694];
                         
                           
prob_bu(:,:,1,1)=prob(1,1)*ones(grid_k_no,ngridkm);
prob_bu(:,:,1,2)=prob(2,1)*ones(grid_k_no,ngridkm);
prob_bu(:,:,2,1)=prob(3,1)*ones(grid_k_no,ngridkm);
prob_bu(:,:,2,2)=prob(4,1)*ones(grid_k_no,ngridkm);

prob_be(:,:,1,1)=prob(1,2)*ones(grid_k_no,ngridkm);
prob_be(:,:,1,2)=prob(2,2)*ones(grid_k_no,ngridkm);
prob_be(:,:,2,1)=prob(3,2)*ones(grid_k_no,ngridkm);
prob_be(:,:,2,2)=prob(4,2)*ones(grid_k_no,ngridkm);

prob_gu(:,:,1,1)=prob(1,3)*ones(grid_k_no,ngridkm);
prob_gu(:,:,1,2)=prob(2,3)*ones(grid_k_no,ngridkm);
prob_gu(:,:,2,1)=prob(3,3)*ones(grid_k_no,ngridkm);
prob_gu(:,:,2,2)=prob(4,3)*ones(grid_k_no,ngridkm);

prob_ge(:,:,1,1)=prob(1,4)*ones(grid_k_no,ngridkm);
prob_ge(:,:,1,2)=prob(2,4)*ones(grid_k_no,ngridkm);
prob_ge(:,:,2,1)=prob(3,4)*ones(grid_k_no,ngridkm);
prob_ge(:,:,2,2)=prob(4,4)*ones(grid_k_no,ngridkm);
