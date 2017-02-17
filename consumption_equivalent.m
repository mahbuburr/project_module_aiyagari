function [ c ] = consumption_equivalent (par, method, sim, U)
%CONSUMPTION EQUIVALENT Calculates the consumption equivalent of a policy
%change
%   par = parameters needed to calculate the consumption equivalent
%   method = specifies whether the simulation method is 'simulation' or 'histogram'
%   sim = matrices with agent's wealth and employment status
%   U = life-time utilities before and and after the policy change
    
    if strcmp(method.sim ,'simulation')
 
        % Calculate the consumption equivalent
        % If value > 1, agents prefer steady state one, if 1 > value > 0, agents prefer
        % policy change ( = steady state two).
        % Consumption equivalent tested against model with policy change:
        c.equivalent = ((U.two.extrap.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.one.extrap.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));

        % Calculate the mean and median consumption equivalent
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
        
        % Sort the consumption equivalent for employed and unemployed in the last period to plot them
        % against wealth
        [c.equivalent_unemp_sorted, sort_index_unemp] = sort(c.equivalent(end,sim.one.e(T,:)==1),'descend');
        [c.equivalent_emp_sorted, sort_index_emp] = sort(c.equivalent(end,sim.one.e(T,:)==2),'descend');
        k_unemp = sim.one.k(T,sim.one.e(T,:)==1);
        k_emp = sim.one.k(T,sim.one.e(T,:)==2);
        
        figure (1)
        plot(k_emp(sort_index_emp),c.equivalent_emp_sorted,'g.',k_unemp(sort_index_unemp),c.equivalent_unemp_sorted,'r.')
        legend('employed','unemployed')
        xlabel('wealth')
        ylabel('consumption equivalent')
        refline (0,1)
       
       
    elseif strcmp(method.sim,'histogram')
        % to be done
    end
end

