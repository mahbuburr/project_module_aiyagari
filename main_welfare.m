clear all
close all
clc

fixed_parameters;
tic
%% Solve the benchmark model
mu_benchmark = 0.5;
sim_e_benchmark = generate_shocks( T, ind_no, L, PI );
k_guess_benchmark = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)
K_lims = [K_rep,grid_k(end)]; % initial limits for bisection method
K_guess_benchmark = (K_lims(1)+K_lims(2))/2;

d1 = 1;
iter = 0;
tic
while d1>1e-6 && iter<50 % loop for aggregate problem
    iter = iter+1;
    
    k_guess_benchmark = solve_individual_problem(mu_benchmark, k_guess_benchmark, K_guess_benchmark);
    [K_demand_benchmark, sim_k_benchmark] = find_dist_agents( sim_e_benchmark, K_guess_benchmark, k_guess_benchmark );
    
    store(1).K_demand(iter) = K_demand_benchmark;
    
    d1 = abs(K_demand_benchmark-K_guess_benchmark)./(1+K_guess_benchmark); % deviation between guess for capital and demanded capital stock

    if K_demand_benchmark>K_guess_benchmark % update limits of bisection interval
        K_lims(1) = K_guess_benchmark;
    else
        K_lims(2) = K_guess_benchmark;
    end
    
    % update guess for aggregate capital stock
    K_guess_benchmark = (K_lims(1)+K_lims(2))/2;
 
    disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K_guess_benchmark),', K_demand: ',num2str(K_demand_benchmark)])
    
end

store(1).mu=mu_benchmark;
store(1).sim_k=sim_k_benchmark(1:50:end,:);
store(1).k_guess = k_guess_benchmark;

%% Solve for the others unemployment benefit levels
period = 5000;
sim_e = generate_shocks( T, ind_no, L, PI, sim_e_benchmark(period,:) );
mu_min = 0.25;
mu_max = 0.90;
mu_n = 10;
grid_mu = linspace(mu_min, mu_max, mu_n);
store(mu_n+1).mu = NaN; %prealocate
count = 2;
parfor nn = 1:mu_n
    mu = grid_mu(nn);
    disp(['mu number: ',num2str(nn)])
    k_guess = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)
    K_lims = [K_rep,grid_k(end)]; % initial limits for bisection method
    K_guess = (K_lims(1)+K_lims(2))/2;
    
    d1 = 1;
    iter = 0;
    tic
    while d1>1e-6 && iter<50 % loop for aggregate problem
        iter = iter+1;
        
        k_guess = solve_individual_problem(mu, k_guess, K_guess);
        [K_demand, sim_k] = find_dist_agents( sim_e, sim_k_benchmark(period,:), k_guess );
        
        store(nn+1).K_demand(iter) = K_demand;
        
        d1 = abs(K_demand-K_guess)./(1+K_guess); % deviation between guess for capital and demanded capital stock
        
        if K_demand>K_guess % update limits of bisection interval
            K_lims(1) = K_guess;
        else
            K_lims(2) = K_guess;
        end
        
        % update guess for aggregate capital stock
        K_guess = (K_lims(1)+K_lims(2))/2;
        
        disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K_guess),', K_demand: ',num2str(K_demand)])
        
    end
    store(nn+1).mu=mu;
    store(nn+1).sim_k=sim_k(1:50:end,:);
    store(nn+1).k_guess = k_guess;
end
store_sim(1).sim_e = sim_e_benchmark;
store_sim(2).sim_e = sim_e;
save('10solutions_riegler2.mat', 'store');
toc




