function [ c, U ] = welfare_effects( par, func, sim, store, K, k, method )
%WELFARE EFFECTS Calculates whether agents prefer a policy change or not
%   Detailed explanation goes here

    if strcmp(method.sim ,'simulation')

        ind_no = size(sim.k,2);
        
        % Get utility levels of all agents at all points in time.
        sim.one.u = func.U(func.C(sim.one.k));
        sim.one.u(func.C(sim.one.k)<0) = -Inf; % Make sure consumption is always positive.
        sim.two.u = func.U(func.C(sim.two.k));
        sim.two.u(func.C(sim.two.k)<0) = -Inf;
        % Calculate the life time utility when simulation has converged (about half
        % the size of the simulation).
        U.lifetime_one = par.beta.^(1:ceil(ind_no/2))*sim.one.u((ceil(ind_no/2)+1):end,:);
        U.lifetime_two = par.beta.^(1:ceil(ind_no/2))*sim.two.u((ceil(ind_no/2)+1):end,:);

        %% Calculate the consumption equivalent
        % If value > 1, agents prefer steady state one, if 1 > value > 0, agents prefer
        % policy change ( = steady state two).
        % Consumption equivalent tested against model with policy change.
        c.equivalent = ((U.lifetime.two*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.lifetime_one*(1-par.sigma)*(1-par.beta)+1)).^(1/(1-par.sigma)); 
        figure(6)
        histogram(c.equivalent)
        xlabel('Consumption equivalents')
        % Calculate average and median of consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
        c.equivalent_mean = mean(c.equivalent);
        c.equivalent_median = median(c.equivalent);

    elseif strcmp(method.sim , 'histogram')
        
        U.lifetime_one = sum(sum(store.one.distribution.*func.U(func.C(k.one.k))));
        U.lifetime_two = sum(sum(store.two.distribution.*func.U(func.C(k.two.k))));
        
        %% Calculate the consumption equivalent
        % If value > 1, agents prefer policy change ( = steady state two), 
        % if 1 > value > 0, agents prefer steady state one.
        % Consumption equivalent tested against model with policy change.
        c.equivalent = ((U.lifetime_two.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.lifetime_one.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma)); 

    end
end