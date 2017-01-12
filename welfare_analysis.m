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
        U.one.guess = func.U(c.one.guess);
        U.one.lifetime = zeros(2,grid.one.k_no); %expected life time utility
        dist=100;
        while dist>1e-8
            EU = par.PI * U.one.lifetime; % Need to rewrite the function because it cannot handle changes in PI
            Unew = U.one.guess' + par.beta * EU;
            dist = max(max(abs(Unew - U.one.lifetime)));
            U.one.lifetime = Unew;
        end
    elseif i==2 % Introduce the policy change here
        par.z = 0.5;
        setup % refresh setup again for new parameter
        [ k.two, c.two, K.two, sim.two, store.two, mat.two, grid.two ] = aiyagari_solver( par, func, method );
        U.two.guess = func.U(c.two.guess);
        U.two.lifetime = zeros(2,grid.two.k_no); %expected life time utility
        dist=100;
        while dist>1e-8
            EU = par.PI * U.two.lifetime; % Need to rewrite the function because it cannot handle changes in PI
            Unew = U.two.guess' + par.beta * EU;
            dist = max(max(abs(Unew - U.two.lifetime)));
            U.two.lifetime = Unew;
        end
    end
end
toc
U.one.lifetime(U.one.lifetime == -Inf) = -999999;
U.two.lifetime(U.two.lifetime == -Inf) = -999999;
ind_no = size(sim.one.k,2);
T = size(sim.one.k,1);
U.one.one = NaN(ceil(T/2),ind_no);
U.two.two = NaN(ceil(T/2),ind_no);
tic
for t = ceil((T+1)/2):T
    U.one.one(t-2500,sim.one.e(t,:)==1) = interp1(grid.one.k, U.one.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
    U.one.one(t-2500,sim.one.e(t,:)==2) = interp1(grid.one.k, U.one.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');

    U.two.two(t-2500,sim.one.e(t,:)==1) = interp1(grid.one.k, U.two.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
    U.two.two(t-2500,sim.one.e(t,:)==2) = interp1(grid.one.k, U.two.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
end
toc
U.one.mean = mean(mean(U.one.one));
U.two.mean = mean(mean(U.two.two));
U.one.median = median(median(U.one.one));
U.two.median = median(median(U.two.two));
% Calculate the consumption equivalent for the two steady states. If c > 1,
% agents prefer the second steady state. If 1 > c > 0, agents prefer the
% first steady state.
[ c ] = consumption_equivalent (par, method, U);
[ k ] = cash_equivalent (par, method, grid, sim, U);
k.equivalent_mean = mean(mean(k.equivalent));
k.equivalent_median = median(median(k.equivalent));