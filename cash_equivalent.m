function [ k ] = cash_equivalent( method, grid, sim, U )
%CASH EQUIVALENT Calculates the capital equivalent of a policy change
%   method = specifies whether the simulation method is 'simulation' or 'histogram'
%   grid = the capital grid to extrapolate the capital compensation of the
%          policy  change
%   sim = matrices with agent's wealth and employment status
%   U = life-time utilities before and and after the policy change

    if strcmp(method.sim ,'simulation')
        ind_no = size(sim.one.k,2);
        T = size(sim.one.k,1);
        k.equivalent = NaN(ceil(T/2),ind_no);
        k.compensated = NaN(ceil(T/2),ind_no);
        
        % Extrapolate the cash equivalent
        for t = 1:ceil(T/2)
            k.compensated(t,sim.one.e(t+ceil(T/2),:)==1) = interp1(U.one.lifetime(1,:), grid.one.k, U.two.extrap(t,sim.one.e(t+ceil(T/2),:)==1), 'linear', 'extrap');
            k.compensated(t,sim.one.e(t+ceil(T/2),:)==2) = interp1(U.one.lifetime(2,:), grid.one.k, U.two.extrap(t,sim.one.e(t+ceil(T/2),:)==2), 'linear', 'extrap');
            k.equivalent(t,:) = k.compensated(t,:) - sim.one.k(t+ceil(T/2),:);   
            k.equivalent_unemployed_mean(t,:) = mean(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
            k.equivalent_unemployed_median(t,:) = median(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
            k.equivalent_employed_mean(t,:) = mean(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
            k.equivalent_employed_median(t,:) = median(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
        end
        
        % Calculate mean and median cash equivalent
        k.equivalent_mean = mean(mean(k.equivalent));
        k.equivalent_median = median(median(k.equivalent));
        
        % Get cash equivalent for employed and unemployed
        k.equivalent_unemployed_mean = mean(k.equivalent_unemployed_mean);
        k.equivalent_unemployed_median = median(k.equivalent_unemployed_median);
        k.equivalent_employed_mean = mean(k.equivalent_employed_mean);
        k.equivalent_employed_median = median(k.equivalent_employed_median);
        
        % Sort the cash equivalent for employed and unemployed in the last period to plot them against
        % wealth
        [k.equivalent_unemp_sorted, sort_index_unemp] = sort(k.equivalent(end,sim.one.e(T,:)==1),'descend');
        [k.equivalent_emp_sorted, sort_index_emp] = sort(k.equivalent(end,sim.one.e(T,:)==2),'descend');
        k_unemp = sim.one.k(T,sim.one.e(T,:)==1);
        k_emp = sim.one.k(T,sim.one.e(T,:)==2);
        output_baseline = 3.3539;
        
        figure (2)
        plot(k_emp(sort_index_emp),k.equivalent_emp_sorted./output_baseline...
            ,'g.',k_unemp(sort_index_unemp),k.equivalent_unemp_sorted./output_baseline,'r.')
        legend('employed','unemployed')
        xlabel('wealth')
        ylabel('cash equivalent / output')
        refline (0,0)
        
        
    elseif strcmp(method.sim, 'histogram')
        % to be done
    end
end

