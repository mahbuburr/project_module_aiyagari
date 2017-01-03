function [ c, U ] = welfare_effects (par, func, method, k, c, K, sim, store, mat, grid)
%WELFARE EFFECTS Summary of this function goes here
%   Detailed explanation goes here
    
    if strcmp(method.sim ,'simulation')
        % Still need to do the transition for the simulation.
        T = size(sim.one.k,1);
        ind_no = size(sim.one.k,2);
        sim.one.c = zeros(T,ind_no); % simulate series of consumption
        sim.two.c = zeros(T,ind_no); % simulate series of consumption
        for t=1:T
            sim.one.c(t,sim.one.e(t,:)==1) = interp1(grid.one.k,c.one.guess(:,1),sim.one.k(t,sim.one.e(t,:)==1),'linear','extrap'); % consumption of currently unemployed
            sim.one.c(t,sim.one.e(t,:)==2) = interp1(grid.one.k,c.one.guess(:,2),sim.one.k(t,sim.one.e(t,:)==2),'linear','extrap'); % consumption of currently employed
            sim.two.c(t,sim.two.e(t,:)==1) = interp1(grid.two.k,c.two.guess(:,1),sim.two.k(t,sim.two.e(t,:)==1),'linear','extrap'); % consumption of currently unemployed
            sim.two.c(t,sim.two.e(t,:)==2) = interp1(grid.two.k,c.two.guess(:,2),sim.two.k(t,sim.two.e(t,:)==2),'linear','extrap'); % consumption of currently employed
        end

        period = ceil(T/2);
        U.one.lifetime = par.beta.^(0:period) * func.U(sim.one.c(period:end,:)); % calculate life-time utility
        U.two.lifetime = par.beta.^(0:period) * func.U(sim.two.c(period:end,:)); % calculate life-time utility
        % Calculate the consumption equivalent
        % If value > 1, agents prefer steady state one, if 1 > value > 0, agents prefer
        % policy change ( = steady state two).
        % Consumption equivalent tested against model with policy change:
        c.equivalent_mean = ((mean(U.two.lifetime).*(1-par.sigma).*(1-par.beta)+1)...
            ./(mean(U.one.lifetime).*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));

        c.equivalent_median = ((median(U.two.lifetime).*(1-par.sigma).*(1-par.beta)+1)...
            ./(median(U.one.lifetime).*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma)); 
        % Calculate average and median of consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
    elseif strcmp(method.sim,'histogram')
        
        for i=1:2 % currently unemployed and employed
            c.one.c(:,i) = max(0,interp1(grid.one.k,c.one.guess(:,i),grid.one.dist,'linear','extrap')); % policy function on grid used for distribution
            c.two.c(:,i) = max(0,interp1(grid.two.k,c.two.guess(:,i),grid.two.dist,'linear','extrap')); % policy function on grid used for distribution
        end
        
        U.one.lifetime = 0; %expected life time utility
        dist = 100;
        while dist>1e-8
            %EU = par.PI * U.one.lifetime;
            Unew = func.U(sum(sum(store.one.distribution .* c.one.c))) + par.beta * U.one.lifetime;
            dist = max(max(abs(Unew - U.one.lifetime)));
            U.one.lifetime = Unew;
        end
        
        U.two.lifetime = 0; %expected life time utility
        dist = 100;
        while dist>1e-8
            %EU = par.PI * U.two.lifetime;
            Unew = func.U(sum(sum(store.one.distribution .* c.two.c))) + par.beta * U.two.lifetime; % multiply with distribution of one to get transition
            dist = max(max(abs(Unew - U.two.lifetime)));
            U.two.lifetime = Unew;
        end
        % Calculate the consumption equivalent
        % If value > 1, agents prefer policy change ( = steady state two), 
        % if 1 > value > 0, agents prefer steady state one.
        % Consumption equivalent tested against model with policy change:
        c.equivalent = ((U.two.lifetime*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.one.lifetime*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));
    end
end

