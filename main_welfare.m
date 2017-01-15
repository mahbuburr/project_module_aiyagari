clear all
close all
clc

fixed_parameters;

%% Solve the benchmark model
mu_benchmark = 0.3;
sim_e_benchmark = generate_shocks( T, ind_no, L, PI );
k_guess = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)
K_guess = (K_lims(1)+K_lims(2))/2;

d1 = 1;
iter = 0;
tic
while d1>1e-6 && iter<50 % loop for aggregate problem
    iter = iter+1;
    
    k_guess = solve_individual_problem(mu_benchmark, k_guess, K_guess);
    [K_demand, sim_k] = find_dist_agents( sim_e_benchmark, K_guess, k_guess );
    
    d1 = abs(K_demand-K_guess)./(1+K_guess); % deviation between guess for capital and demanded capital stock
    
    store.K_guess(iter)=K_guess;
    store.K_next(iter)=K_demand;
    
    if K_demand>K_guess % update limits of bisection interval
        K_lims(1) = K_guess;
    else
        K_lims(2) = K_guess;
    end
    
    % update guess for aggregate capital stock
    K_guess = (K_lims(1)+K_lims(2))/2;
 
    disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K_guess),', K_demand: ',num2str(K_demand)])
    
end


