function [ c ] = consumption_equivalent (par, method, U, sim)
%WELFARE EFFECTS Summary of this function goes here
%   Detailed explanation goes here
    
    if strcmp(method.sim ,'simulation')
 
        % Calculate the consumption equivalent
        % If value > 1, agents prefer steady state one, if 1 > value > 0, agents prefer
        % policy change ( = steady state two).
        % Consumption equivalent tested against model with policy change:
        c.equivalent = ((U.two.extrap.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.one.extrap.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));

        % Calculate average and median of consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
        c.equivalent_mean = mean(mean(c.equivalent));
        c.equivalent_median = median(median(c.equivalent));
        
        % Get consumption equivalent for employed and unemployed
        T = size(sim.one.e,1);
        for t=1:ceil(T/2)
            c.equivalent_unemployed_mean(t,:) = mean(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
            c.equivalent_unemployed_median(t,:) = median(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
            c.equivalent_employed_mean(t,:) = mean(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
            c.equivalent_employed_median(t,:) = median(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
        end
        c.equivalent_unemployed_mean = mean(c.equivalent_unemployed_mean);
        c.equivalent_unemployed_median = median(c.equivalent_unemployed_median);
        c.equivalent_employed_mean = mean(c.equivalent_employed_mean);
        c.equivalent_employed_median = median(c.equivalent_employed_median);
        
        
%         [c.equivalent_unemp_sorted, sort_index_unemp] = sort(c.equivalent(end,sim.one.e(T,:)==1),'descend');
%         [c.equivalent_emp_sorted, sort_index_emp] = sort(c.equivalent(end,sim.one.e(T,:)==2),'descend');
%         k_unemp = sim.one.k(T,sim.one.e(T,:)==1);
%         k_emp = sim.one.k(T,sim.one.e(T,:)==2);
%         figure (1)
%         plot(k_unemp(sort_index_unemp),c.equivalent_unemp_sorted,'g',k_emp(sort_index_emp),c.equivalent_emp_sorted,'r')
%         legend('unemployed','employed')
%         xlabel('wealth')
%         ylabel('consumption equivalent')
%         refline (0,1)
        
    elseif strcmp(method.sim,'histogram')
        % to be done
        
%         for i=1:2 % currently unemployed and employed
%             c.one.c(:,i) = max(0,interp1(grid.one.k,c.one.guess(:,i),grid.one.dist,'linear','extrap')); % consumption policy function on grid used for distribution
%             c.two.c(:,i) = max(0,interp1(grid.two.k,c.two.guess(:,i),grid.two.dist,'linear','extrap')); % policy function on grid used for distribution
%         end
%         
%         U.one.lifetime = 0; %expected life time utility
%         dist = 100;
%         while dist>1e-8
%             %EU = par.PI * U.one.lifetime;
%             Unew = func.U(sum(sum(store.one.distribution .* c.one.c))) + par.beta * U.one.lifetime;
%             dist = max(max(abs(Unew - U.one.lifetime)));
%             U.one.lifetime = Unew;
%         end
%         
%         U.two.lifetime = 0; %expected life time utility
%         dist = 100;
%         while dist>1e-8
%             %EU = par.PI * U.two.lifetime;
%             Unew = func.U(sum(sum(store.two.distribution .* c.two.c))) + par.beta * U.two.lifetime; % multiply with distribution of one to get transition
%             dist = max(max(abs(Unew - U.two.lifetime)));
%             U.two.lifetime = Unew;
%         end
%         % Calculate the consumption equivalent
%         % If value > 1, agents prefer policy change ( = steady state two), 
%         % if 1 > value > 0, agents prefer steady state one.
%         % Consumption equivalent tested against model with policy change:
%         c.equivalent = ((U.two.lifetime*(1-par.sigma).*(1-par.beta)+1)...
%             ./(U.one.lifetime*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));
    end
end

