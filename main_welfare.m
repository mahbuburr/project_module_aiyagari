clear all
close all
clc

fixed_parameters;

%% Solve the benchmark model
mu_benchmark = 0.3;
sim_e_benchmark = generate_shocks( T, ind_no, L, PI );
k_guess = solve_individual_problem( mu, L, grid_k_no, grid_k, k_guess, K_guess, PI, beta, delta, mat_k,  k_min);

% %% solve for general equilibrium
% HH_method = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
% sim_method = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
% agg_method = 'bisection'; % 'bisection' to use bisection method, gradual updating otherwise
% 


% initial guesses


