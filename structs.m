% s = struct('mu',mu_benchmark,'k_sim',k_guess_benchmark);
% b = struct('mu',mu_benchmark,'k_sim',k_guess_benchmark);
% 
% field = 's';
% value = {mu_benchmark;
%          k_guess_benchmark};
% s = struct(field,value)

store(1).mu=mu_benchmark;
store(1).sim_k=sim_k_benchmark;
store(1).k_guess = k_guess_benchmark;
store(1).K_demand = K_demand_benchmark;

store(2).mu=mu;
store(2).sim_k=sim_k;
store(2).k_guess = k_guess;
store(2).K_demand = K_demand;

% s(1) = struct('mu',mu_benchmark, 'sim_k',sim_k_benchmark);