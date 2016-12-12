clear
close all

parameters % load paramaters
method.sim = 'simulation'
tic
[k, c, K, sim, store]= aiyagari_solver(par, grid, K, k, func, method, mat);
toc
