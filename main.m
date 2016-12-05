clear
close all

parameters % load paramaters

tic % Set timer
% Solve the Aiyagariy model by fixed-point iteration.
% Function in- and output are structures.
[k, c, K, sim, store]= aiyagari_solver(par, grid, K, k, func, method, mat);
toc