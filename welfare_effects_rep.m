function [ c, U ] = welfare_effects_rep (par, func, method, k, K, sim, store, mat, grid)
%WELFARE EFFECTS REP Calculates whether agents prefer a frictionless, representative agent model or not
%   Detailed explanation goes here

    c.guess = (1+func.r(K.guess)-par.delta)*mat.k-k.guess+mat.income(K.guess); % Get consumption policy function
    c.guess = max(0, c.guess);
    
    if strcmp(method.sim ,'simulation')

        
        T = size(sim.k);
        ind_no = size(sim.k,2);
        sim.c = zeros(T,ind_no); % simulate series of consumption
        for t=1:T
            sim.c(t,sim.e(t,:)==1) = interp1(grid.k,c.guess(:,1),sim.k(t,sim.e(t,:)==1),'linear','extrap'); % consumption of currently unemployed
            sim.c(t,sim.e(t,:)==2) = interp1(grid.k,c.guess(:,2),sim.k(t,sim.e(t,:)==2),'linear','extrap'); % consumption of currently employed
        end

        period = ceil(T/2);
        U.lifetime = par.beta.^(0:period) * func.U(sim.c(period:end,:)); % calculate life-time utility
        U.lifetime_rep = sum(par.beta.^(0:period) * func.U(func.C(K.rep)));
       
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
%         a = sim.c(period:end,:);
%         figure(7)
%         histogram(a(:),100)
        

    elseif strcmp(method.sim , 'histogram')
   
        for i=1:2 % currently unemployed and employed
            c.c(:,i) = max(0,interp1(grid.k,c.guess(:,i),grid.dist,'linear','extrap')); % policy function on grid used for distribution
        end
        
        U.lifetime = 0; %expected life time utility
        dist = 100;
        while dist>1e-8
            EU = par.PI * U.lifetime;
            Unew = func.U(sum(sum(store.distribution .* c.c))) + par.beta * EU;
            dist = max(max(abs(Unew - U.lifetime)));
            U.lifetime = Unew;
        end
        
        
        U.lifetime_rep = 0; %expected life time utility representative agent
        dist = 100;
        while dist>1e-8
            EU = par.PI * U.lifetime_rep;
            Unew = func.U(sum(sum(store.distribution .* (func.C(K.rep))))) + par.beta * EU;
            dist = max(max(abs(Unew - U.lifetime_rep)));
            U.lifetime_rep = Unew;
        end

        %% Calculate the consumption equivalent
        % If value > 1, agents prefer standard model, if 1 > value > 0, agents prefer
        % Aiyagari model.
        % Consumption equivalent tested against frictionless benchmark model.
        c.equivalent_bench = ((U.lifetime_rep.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.lifetime.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma)); 
    end
end

