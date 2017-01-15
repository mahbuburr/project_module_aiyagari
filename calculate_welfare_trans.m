clear all
close all
load('10solutions_riegler_trans.mat');
period = 5000; 
% Calculate expected life time utility 
% for the benchmark economy in steady state
U_benchmark = expec_lifetime_utility( store(1).k_guess, store(1).K_demand, store(1).mu, shocks(1).sim_e, store(1).sim_k, period );

%% Bla
time = 1:50:T;
U = NaN(length(time), ind_no);
U(1, :) = U_benchmark;

% Calculate expected life time utility 
% for the comparison economy in transition

parfor t = 1:length(time)
    Uaux = expec_lifetime_utility( k_guess, K_demand, mu, sim_e, sim_k, time(t) );
    U(t+1,:) = Uaux;
end