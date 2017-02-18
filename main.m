clear
close all

parameters % load paramaters
% par.mu = 0.01;
setup % load setup
% method.HH = 'FPend'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
% method.agg = 'bisection'; % 'bisection' to use bisection method, gradual updating otherwise
tic
[k, c, K, sim, store]= aiyagari_solver(par, func, method);
toc
