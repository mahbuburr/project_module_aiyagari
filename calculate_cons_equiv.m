function [ cons_equiv, cons_equiv_mean, cons_equiv_median ] = calculate_cons_equiv( mu_n, benchmark, store, name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load(strcat('fixed_parameters_',name,'.mat'));

cons_equiv = ((store(benchmark).U*(1-sigma)*(1-betta)+1)./(store(mu_n).U*(1-sigma)*(1-betta)+1)).^(1/(1-sigma));
    
cons_equiv_mean = mean(cons_equiv);
cons_equiv_median = median(cons_equiv);
end

