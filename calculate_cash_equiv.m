function [ cash_equiv, cash_equiv_agg] = calculate_cash_equiv( mu_n, benchmark, store, name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load(strcat('fixed_parameters_',name,'.mat'));

%ASK HERE
cash_equiv= NaN(1,ind_no);
for ind = 1:ind_no
    cash_equiv(ind) = interp1(store(mu_n).U, grid_k, store(benchmark).U) - store(benchmark).sim_k;
end

cash_equiv_agg = sum(cash_equiv);
end