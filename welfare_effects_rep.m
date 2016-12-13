function [ c, U ] = welfare_effects_rep( par, func, sim, store, K, k, method )
%WELFARE EFFECTS REP Calculates whether agents prefer a frictionless, representative agent model or not
%   Detailed explanation goes here

    if strcmp(method.sim ,'simulation')

        ind_no = size(sim.k,2);
     
        % Calculate the life time utility when simulation has converged (about half
        % the size of the simulation).
        U.lifetime = par.beta.^(1:ceil(ind_no/2))*func.U(func.C(sim.k((ceil(ind_no/2)+1):end,:)));
        U.lifetime_rep = sum(par.beta.^(1:ceil(ind_no/2))*func.U(func.C(K.rep)));

        %% Calculate the consumption equivalent
        % If value > 1, agents prefer standard model, if 1 > value > 0, agents prefer
        % Aiyagari model.
        % Consumption equivalent tested against frictionless benchmark model.
        c.equivalent_bench = ((U.lifetime_rep.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.lifetime.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma)); 
        figure(6)
        histogram(c.equivalent_bench)
        xlabel('Consumption equivalents')
        % Calculate average and median of consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
        c.equivalent_bench_mean = mean(c.equivalent_bench);
        c.equivalent_bench_median = median(c.equivalent_bench);
        
        % Just to see how consumption looks overall
        a = func.C(sim.k((ceil(ind_no/2)+1):end,:));
        figure(7)
        histogram(a(:),100)
        

    elseif strcmp(method.sim , 'histogram')
        
        U.lifetime = func.U(sum(sum(store.distribution.*(func.C(k.k))))); % Aggregate utility
        U.lifetime_rep = func.U(func.C(K.rep));

        %% Calculate the consumption equivalent
        % If value > 1, agents prefer standard model, if 1 > value > 0, agents prefer
        % Aiyagari model.
        % Consumption equivalent tested against frictionless benchmark model.
        c.equivalent_bench = ((U.lifetime_rep.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.lifetime.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma)); 
    end
end

