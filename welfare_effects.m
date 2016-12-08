function [ c ] = welfare_effects( par, func, sim, store, K, k, method )
%WELFARE EFFECTS Calculates whether agents prefer a policy change or not
%   Detailed explanation goes here

    % Define function to calculate utility of consumption.
    func.U = @(c) (c.^(1-par.sigma))./(1-par.sigma);

    if strcmp(method.sim ,'simulation')

        ind_no = size(sim.k,2);
        
        % Get utility levels of all agents at all points in time.
        sim.u = func.U(func.C(sim.k));
        sim.u(func.C(sim.k)<0) = -Inf; % Make sure consumtion is always positive.
        % Calculate the life time utility when simulation has converged (about half
        % the size of the simulation).
        lifetime_U = par.beta.^(1:ceil(ind_no/2))*sim.u((ceil(ind_no/2)+1):end,:);
        lifetime_U_rep = sum(par.beta.^(1:ceil(ind_no/2))*func.U(func.C(K.rep)));

        %% Calculate the consumption equivalent
        % If value>1, agents prefer old steady state, if 1>value>0, agents prefer
        % policy change.
        % Consumption equivalent tested against frictionless benchmark model.
        c.equivalent_bench = ((lifetime_U.*(1-par.sigma).*(1-par.beta)+1)...
            ./(lifetime_U_rep*(1-par.sigma)*(1-par.beta)+1)).^(1/(1-par.sigma)); 
        figure(6)
        histogram(c.equivalent_bench)
        xlabel('Consumption equivalents')
        % Calculate average and median of consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
        c.equivalent_bench_mean = mean(c.equivalent_bench);
        c.equivalent_bench_median = median(c.equivalent_bench);

    elseif strcmp(method.sim , 'histogram')
        
        lifetime_U = store.distribution'*func.U(func.C(k.k));
        lifetime_U_rep = func.U(func.C(K.rep));

        %% Calculate the consumption equivalent
        % If value>1, agents prefer old steady state, if 1>value>0, agents prefer
        % policy change.
        % Consumption equivalent tested against frictionless benchmark model.
        c.equivalent_bench = ((diag(lifetime_U).*(1-par.sigma).*(1-par.beta)+1)...
            ./(lifetime_U_rep*(1-par.sigma)*(1-par.beta)+1)).^(1/(1-par.sigma)); 

        % Calculate the average of the consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
        c.equivalent_bench_mean = (1-sim.L)*c.equivalent_bench(1) + sim.L*c.equivalent_bench(2);

    end
end

