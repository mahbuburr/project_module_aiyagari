function [ k ] = cash_equivalent( par, method, grid, sim, U )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(method.sim ,'simulation')
        ind_no = size(sim.one.k,2);
        T = size(sim.one.k,1);
        k.equivalent = NaN(ceil(T/2),ind_no);
        k.compensated = NaN(ceil(T/2),ind_no);
        
        for t = 1:ceil(T/2)
            k.compensated(t,sim.one.e(t+ceil(T/2),:)==1) = interp1(U.one.lifetime(1,:), grid.one.k, U.two.extrap(t,sim.one.e(t+ceil(T/2),:)==1), 'linear', 'extrap');
            k.compensated(t,sim.one.e(t+ceil(T/2),:)==2) = interp1(U.one.lifetime(2,:), grid.one.k, U.two.extrap(t,sim.one.e(t+ceil(T/2),:)==2), 'linear', 'extrap');
            k.equivalent(t,:) = k.compensated(t,:) - sim.one.k(t+ceil(T/2),:);   
            k.equivalent_unemployed_mean(t,:) = mean(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
            k.equivalent_unemployed_median(t,:) = median(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
            k.equivalent_employed_mean(t,:) = mean(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
            k.equivalent_employed_median(t,:) = median(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
        end
        
        k.equivalent_mean = mean(mean(k.equivalent));
        k.equivalent_median = median(median(k.equivalent));
        
        % Get cash equivalent for employed and unemployed
        k.equivalent_unemployed_mean = mean(k.equivalent_unemployed_mean);
        k.equivalent_unemployed_median = median(k.equivalent_unemployed_median);
        k.equivalent_employed_mean = mean(k.equivalent_employed_mean);
        k.equivalent_employed_median = median(k.equivalent_employed_median);
        
    elseif strcmp(method.sim, 'histogram')
        % to be done
end

