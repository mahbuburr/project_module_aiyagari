function [ U ] = calculate_welfare_trans( mu_n, store, shocks )
% Calculate expected life time utility 
% for the benchmark economy in steady state
fixed_parameters;
time = size(store(1).sim_k,1)+1;
U_benchmark = expec_lifetime_utility_trans( store(1).k_guess, store(1).K_demand, store(1).mu, shocks(1).sim_e, store(1).sim_k, 1 );

U = NaN(time, ind_no);
U(1, :) = U_benchmark;

% Calculate expected life time utility 
% for the comparison economy in transition
for t = 1:time-1
    Uaux = expec_lifetime_utility_trans( store(mu_n).k_guess, store(mu_n).K_demand, store(mu_n).mu, shocks(2).sim_e, store(mu_n).sim_k, t);
    U(t+1,:) = Uaux;
end

end