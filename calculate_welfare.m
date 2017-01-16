function [ U, U_mean ] = calculate_welfare( mu_n, store, shocks, name )
% Calculate expected life time utility 
% for the benchmark economy in steady state
load(strcat('fixed_parameters_',name,'.mat'));

e = shocks.sim_e(5000,:); 
elu = expec_lifetime_utility( store(mu_n).k_guess, store(mu_n).K_demand, store(mu_n).mu, name);
U = NaN(1,5000);
U(e==1) = interp1(grid_k,elu(1,:),store(mu_n).sim_k(e==1),'linear','extrap');
U(e==2) = interp1(grid_k,elu(2,:),store(mu_n).sim_k(e==2),'linear','extrap');
U_mean = mean(U);
end