function [ cons_equiv, cons_equiv_mean, cons_equiv_median ] = calculate_cons_equiv( mu_n, benchmark, store )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fixed_parameters;

cons_equiv = ((store(mu_n).U*(1-sigma)*(1-betta)+1)./(store(benchmark).U*(1-sigma)*(1-betta)+1)).^(1/(1-sigma));
    
cons_equiv_mean = mean(cons_equiv);
cons_equiv_median = median(cons_equiv);
end

