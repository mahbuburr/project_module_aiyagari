clear
close all

parameters % load paramaters
setup % load setup

% Change simulation method to histogram if you want.
%method.sim = 'histogram';

%% Welfare effects of two Aiyagari models with different parameters.

% Default parameterization is set by the calibration of Markus Riegler. If
% you change parameters type "par.<parameter> =" and then "setup" before
% you call the function. If you want to compare one steady state with the
% default one, change parameters for i=2 and leave i=1 as the default one.
% Do not change beta or sigma!
mu_min = 0.05;
mu_max = 0.8;
mu_n = 20;
mu = linspace(mu_min, mu_max, mu_n);
tic
for i=1:2 % Get steady state utility for the model with two different parameters
    if i==1 
        [ k.one, c.one, K.one, sim.one, store.one, mat.one, grid.one ] = aiyagari_solver( par, func, method );
        U.one.guess = func.U(c.one.guess);
        U.one.lifetime = zeros(2,grid.one.k_no); %expected life time utility
        dist=100;
        while dist>1e-8
            EU = par.PI * U.one.lifetime; 
            Unew = U.one.guess' + par.beta * EU;
            dist = max(max(abs(Unew - U.one.lifetime)));
            U.one.lifetime = Unew;
        end
        U.one.lifetime(U.one.lifetime == -Inf) = -999999; % Get rid of -Inf for negative consumption for extrapolation
        ind_no = size(sim.one.k,2);
        T = size(sim.one.k,1);
        U.one.extrap = NaN(ceil(T/2),ind_no);
        for t = ceil((T+1)/2):T % Extrapolated life time utility with transition
            U.one.extrap(t-ceil(T/2),sim.one.e(t,:)==1) = interp1(grid.one.k, U.one.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
            U.one.extrap(t-ceil(T/2),sim.one.e(t,:)==2) = interp1(grid.one.k, U.one.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
        end
    elseif i==2 
        for ii=5:size(mu,2)
            par.mu = mu(ii)
            method.HH = 'FP'; 
            method.sim = 'simulation'; 
            method.agg = 'bisection';
            setup % refresh setup for new parameter
            [ k.two, c.two, K.two, sim.two, store.two, mat.two, grid.two ] = aiyagari_solver( par, func, method );
            U.two.guess = func.U(c.two.guess);
            U.two.lifetime = zeros(2,grid.two.k_no); %expected life time utility
            dist=100;
            while dist>1e-8
                EU = par.PI * U.two.lifetime; 
                Unew = U.two.guess' + par.beta * EU;
                dist = max(max(abs(Unew - U.two.lifetime)));
                U.two.lifetime = Unew;
            end
            U.two.lifetime(U.two.lifetime == -Inf) = -999999; % Get rid of -Inf for extrapolation
            ind_no = size(sim.one.k,2);
            T = size(sim.one.k,1);
            U.two.extrap = NaN(ceil(T/2),ind_no);
            for t = ceil((T+1)/2):T % Extrapolated life time utility with transition
                U.two.extrap(t-ceil(T/2),sim.one.e(t,:)==1) = interp1(grid.one.k, U.two.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
                U.two.extrap(t-ceil(T/2),sim.one.e(t,:)==2) = interp1(grid.one.k, U.two.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
            end
            % Calculate the consumption equivalent for the two steady states. If c > 1,
            % agents prefer the second steady state. If 1 > c > 0, agents prefer the
            % first steady state.
            [ c ] = consumption_equivalent (par, method, U, sim);

            % Calculate the cash equivalent for the two steady states.
            [ k ] = cash_equivalent (par, method, grid, sim, U);

            keep.c.equivalent_mean = c.equivalent_mean;
            keep.c.equivalent_median = c.equivalent_median;
            keep.c.equivalent_unemployed_mean = c.equivalent_unemployed_mean;
            keep.c.equivalent_unemployed_median = c.equivalent_unemployed_median;
            keep.c.equivalent_employed_mean = c.equivalent_employed_mean;
            keep.c.equivalent_employed_median = c.equivalent_employed_median;
            keep.k.equivalent_mean = k.equivalent_mean;
            keep.k.equivalent_median = k.equivalent_median;
            keep.k.equivalent_unemployed_mean = k.equivalent_unemployed_mean;
            keep.k.equivalent_unemployed_median = k.equivalent_unemployed_median;
            keep.k.equivalent_employed_mean = k.equivalent_employed_mean;
            keep.k.equivalent_employed_median = k.equivalent_employed_median;
            keep.K = K.two.guess;
            filename = ['baseline_mu_' num2str(ii) '.mat'];
            save (filename, '-struct','keep');
        end
    end
end
toc


