period = 5000; 
% Calculate expected life time utility 
% for the benchmark economy in steady state
U_benchmark = expec_lifetime_utility( k_guess_benchmark, K_guess_benchmark, mu_benchmark, sim_e_benchmark, sim_k_benchmark, period );

%% Bla
time = 1:50:T;
U = NaN(length(time), ind_no);
U(1, :) = U_benchmark;

% Calculate expected life time utility 
% for the comparison economy in transition

parfor t = 1:length(time)
    Uaux = expec_lifetime_utility( k_guess, K_guess, mu, sim_e, sim_k, time(t) );
    U(t+1,:) = Uaux;
end