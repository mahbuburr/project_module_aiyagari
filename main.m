% clear
% close all

parameters % load paramaters
setup % load setup
method.sim = 'simulation';
tic
[k, c, K, sim, store]= aiyagari_solver(par, func, method);
toc
