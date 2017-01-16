clear
close all
clc

fixed_parameters;
tic

%% Solve for the unemployment benefit levels
period = 5000;
sim_e = generate_shocks( T, ind_no, L, PI);
mu_min = 0.30;
mu_max = 0.99;
mu_n = 20;
grid_mu = linspace(mu_min, mu_max, mu_n);
store(mu_n).mu = NaN; %prealocate
for nn = 1:mu_n
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
        [K_demand, sim_k] = find_dist_agents( sim_e, K_guess, k_guess );
        
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
    store(nn).mu=mu;
    store(nn).sim_k=sim_k(period,:);
    store(nn).k_guess = k_guess;
    store(nn).K_demand = K_demand;
    store(nn).K_demand = K_guess;
end
shocks.sim_e = sim_e;
save('20solutions_riegler.mat', 'store', 'shocks');
toc




