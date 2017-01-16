function [ U ] = calculate_welfare_trans( mu_n, store, shocks )
% Calculate expected life time utility 
% for the benchmark economy in steady state
fixed_parameters;
time = size(store(1).sim_k,1)+1;
e = shocks(1).sim_e_benchmark(5000,:); 
expec_lifetime_utility = expec_lifetime_utility_trans( store(1).k_guess, store(1).K_demand, store(1).mu);
U_benchmark = NaN(1,ind_no);
U_benchmark(e(individual)==1) = interp1(grid_k,expec_lifetime_utility(e(individual)==1,:),sim_k(1,e(individual)==1),'linear','extrap');
U_benchmark(e(individual)==2) = interp1(grid_k,expec_lifetime_utility(e(individual)==2,:),sim_k(1,e(individual)==2),'linear','extrap');


U = NaN(time, ind_no);
U(1, :) = U_benchmark;

% Calculate expected life time utility 
% for the comparison economy in transition
for t = 1:time-1
    expec_lifetime_utility = expec_lifetime_utility_trans( store(mu_n).k_guess, store(mu_n).K_demand, store(mu_n).mu);
    e = shocks(2).sim_e_benchmark(t,:); 
    U_aux(e(individual)==1) = interp1(grid_k,expec_lifetime_utility(e(individual)==1,:),sim_k(t,e(individual)==1),'linear','extrap');
    U_aux(e(individual)==2) = interp1(grid_k,expec_lifetime_utility(e(individual)==2,:),sim_k(t,e(individual)==2),'linear','extrap');
    U(t+1,:) = Uaux;
end

end