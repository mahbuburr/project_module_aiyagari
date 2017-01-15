clear all
close all
clc

fixed_parameters;
sim_e = generate_shocks( T, ind_no, L, PI );

mu_benchmark = 0.5;

%% solve for general equilibrium
HH_method = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
sim_method = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
agg_method = 'bisection'; % 'bisection' to use bisection method, gradual updating otherwise

K_lims = [K_rep,grid_k(end)]; % initial limits for bisection method

% initial guesses
k_guess = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)

K_guess = (K_lims(1)+K_lims(2))/2;
