% s = struct('mu',mu_benchmark,'k_sim',k_guess_benchmark);
% b = struct('mu',mu_benchmark,'k_sim',k_guess_benchmark);
% 
% field = 's';
% value = {mu_benchmark;
%          k_guess_benchmark};
% s = struct(field,value)


sto(1).mu=mu_benchmark;
sto(1).sim_k=sim_k_benchmark;
sto(1).k_guess = k_guess_benchmark;
sto(1).K_demand = K_demand_benchmark;
% save('proba.mat', 'sto');
% clear sto

sto(2).mu=mu;
sto(2).sim_k=sim_k;
sto(2).k_guess = k_guess;
sto(2).K_demand = K_demand;
% save('proba.mat', 'sto', '-append');
% clear sto
% s(1) = struct('mu',mu_benchmark, 'sim_k',sim_k_benchmark);