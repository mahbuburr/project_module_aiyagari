clear
close all

parameters % load paramaters
setup % load setup
method.HH = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
method.agg = 'bisectio'; % 'bisection' to use bisection method, gradual updating otherwise
tic
[k, c, K, sim, store]= aiyagari_solver(par, func, method);
toc
