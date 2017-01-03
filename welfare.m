clear
close all

parameters % load paramaters
setup % load setup

% Change simulation method to simulation if you want, it takes much longer
% though.
method.sim = 'simulation';

%% Welfare analysis against frictionless representative agent model

% tic % Set timer
% % Solve the Aiyagariy model by fixed-point iteration.
% % Function in- and output are structures.
% [k, c, K, sim, store, mat, grid]= aiyagari_solver(par, func, method);
% toc
% 
% % Compare welfare against frictionless model
% [c, U] = welfare_effects_rep(par, func, method, k, c, K, sim, store, mat, grid);

%% Welfare effects of two Aiyagari models with different parameters.

% Default parameterization is set by the calibration of Markus Riegler. If
% you change parameters type "par.<parameter> =" and then "setup" before
% you call the function. If you want to compare one steady state with the
% default one, change parameters for i=2 and leave i=1 as the default one.
% Do not change beta or sigma!

tic
for i=1:2 % Get steady state utility for the model with two different parameters
    if i==1 
        [ k.one, c.one, K.one, sim.one, store.one, mat.one, grid.one ] = aiyagari_solver( par, func, method );
    elseif i==2
        par.z = 0.5;
        setup % refresh setup again for new parameter
        [ k.two, c.two, K.two, sim.two, store.two, mat.two, grid.two ] = aiyagari_solver( par, func, method );
    end
end
toc

% Calculate the consumption equivalent for the two steady states. If c > 1,
% agents prefer the second steady state. If 1 > c > 0, agents prefer the
% first steady state.
[ c, U ] = welfare_effects (par, func, method, k, c, K, sim, store, mat, grid);
