%% Welfare analysis of two Aiyagari models with different parameters
clear
close all

parameters % load paramaters
setup % load setup

% Create grid of unemployment benefits for the analysis
mu_min = 0.01;
mu_max = 0.6;
mu_n = 20;
mu = linspace(mu_min, mu_max, mu_n);
mu(7) = mu(6); % the original point does not converge
mu(6) = 0.15; % the original point does not converge
mu(8) = 0.18; % the original point does not converge

% Corresponding PI_UE grid for analysis 
PI_UE_grid = [0.418006431, 0.413867753, 0.409823146, 0.405844156, 0.401954115,0.4, 0.398125746, 0.396341463, 0.38708909, 0.383537395, 0.380061395, 0.376636922, 0.373284328, 0.369980363, 0.366744717, 0.363555009, 0.360430298, 0.35734902, 0.354329636, 0.351351351];
% Default parameterization is set by the calibration of Krussel and Smith (1998) for a US recession.
% If you change parameters, type "par.<parameter> =" and then "setup" before
% you call the function. Do not change beta or sigma!
% If you want to compare one steady state with the baseline, change
% parameters for i=2 and leave i=1 as the baseline.

% Get steady state utility and welfare measures for the model with two different parameters
for i=1:2
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
        for t = ceil((T+1)/2):T % Extrapolate life time utility
            U.one.extrap(t-ceil(T/2),sim.one.e(t,:)==1) = interp1(grid.one.k, U.one.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
            U.one.extrap(t-ceil(T/2),sim.one.e(t,:)==2) = interp1(grid.one.k, U.one.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
        end
    elseif i==2 
        for ii=1%:size(mu,2) 
            tic
            iteration = ii
            par.mu = mu(ii);
            par.PI_UE = PI_UE_grid(ii);
            %method.HH = 'FPend'; % Depending on the mu, you might have to change it to 'FP' or 'FPend' to converge
            %method.agg = 'bisection'; % Depending on the mu, you might have to change it to 'bisection' or 'bisectio' to converge
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
            for t = ceil((T+1)/2):T % Extrapolate life time utility with idiosyncratic transition
                U.two.extrap(t-ceil(T/2),sim.one.e(t,:)==1) = interp1(grid.one.k, U.two.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
                U.two.extrap(t-ceil(T/2),sim.one.e(t,:)==2) = interp1(grid.one.k, U.two.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
            end
            % Calculate the consumption equivalentof the policy reform. If c.equivalent > 1,
            % agents prefer the policy change. If 1 > c.equivalent > 0, agents prefer the
            % baseline model.
            [ c ] = consumption_equivalent (par, method, sim, U);

            % Calculate the cash equivalent of the policy change. If
            % k.equivalent > 0, agents prefer the policy change. If
            % k.equivalent < 0, agents prefer the baseline model.
            [ k ] = cash_equivalent (method, grid, sim, U);

%             % Construct a variable of values to be stored
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
            filename = ['adapting_transitions_2_mu_' num2str(ii) '.mat'];
            %filename = ['adapting_transitions_mu_' num2str(ii) '.mat']; %Changing transitions
            %filename = ['baseline_mu_' num2str(ii) '.mat']; %Baseline
            save (filename, '-struct','keep'); % Save the values
            toc
        end
    end
end