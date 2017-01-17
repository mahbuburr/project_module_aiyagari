function [ cash_equiv, cash_equiv_agg] = calculate_cash_equiv( mu_n, benchmark, store, sim_e, name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load(strcat('fixed_parameters_',name,'.mat'));


cash_equiv= NaN(1,ind_no);
elu_bench = expec_lifetime_utility( store(benchmark).k_guess, store(benchmark).K_demand, store(benchmark).mu, name);
elu_change = expec_lifetime_utility( store(mu_n).k_guess, store(mu_n).K_demand, store(mu_n).mu, name);
for ind = 1:ind_no
    U = interp1(grid_k,elu_bench(sim_e(ind),:),store(benchmark).sim_k(ind),'linear','extrap');
    cash_equiv(ind) = interp1(elu_change(sim_e(ind),:), grid_k, U) - store(benchmark).sim_k(ind);
end

cash_equiv_agg = sum(cash_equiv);
end