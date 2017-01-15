function [ sim_e ] = generate_shocks( T, ind_no, L, PI, sim_e_initial )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
rng('default') % reset random number generator

switch nargin
    case 4
        sim_e = ones(T,ind_no); % simulated employment status
        sim_e(1,1:round(L*ind_no))=2; % initial individuals that are employed
    case 5
        sim_e = NaN(T, ind_no);
        sim_e(1,:) = sim_e_initial;
end

sim_shock = rand(T,ind_no); % shocks for employment transition

for t=2:T
    sim_e(t,sim_e(t-1,:)==1) = 1+(sim_shock(t,sim_e(t-1,:)==1)<=PI(1,2)); % new employment status of previously unemployed
    sim_e(t,sim_e(t-1,:)==2) = 1+(sim_shock(t,sim_e(t-1,:)==2)<=PI(2,2)); % new employment status of previously employed
end

end

