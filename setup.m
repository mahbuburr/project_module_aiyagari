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