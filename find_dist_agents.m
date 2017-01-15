function [ K_demand, sim_k ] = find_dist_agents( sim_e, k_initial, k_guess )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
fixed_parameters;

sim_k = NaN(T,ind_no); % simulated values of capital stock

sim_k(1,:) = k_initial; % initial capital holdings

for t=2:T
    sim_k(t,sim_e(t,:)==1) = interp1(grid_k,k_guess(:,1),sim_k(t-1,sim_e(t,:)==1),'linear','extrap'); % capital demand of currently unemployed
    sim_k(t,sim_e(t,:)==2) = interp1(grid_k,k_guess(:,2),sim_k(t-1,sim_e(t,:)==2),'linear','extrap'); % capital demand of currently employed
end

K_demand = mean(mean(sim_k(ceil(T/2):end,:))); % average capital holdings over second half of sample

end

