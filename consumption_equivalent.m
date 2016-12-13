function [ c ] = consumption_equivalent( U, par, method )
%WELFARE EFFECTS Summary of this function goes here
%   Detailed explanation goes here
    
    if strcmp(method.sim ,'simulation')
        % Calculate the consumption equivalent
        % If value > 1, agents prefer steady state one, if 1 > value > 0, agents prefer
        % policy change ( = steady state two).
        % Consumption equivalent tested against model with policy change:
        c.equivalent_mean = ((mean(U.two).*(1-par.sigma).*(1-par.beta)+1)...
            ./(mean(U.one).*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));

        c.equivalent_median = ((median(U.two).*(1-par.sigma).*(1-par.beta)+1)...
            ./(median(U.one).*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma)); 
        % Calculate average and median of consumption equivalent. If this aggregate
        % >1, there exists a (lump-sum) redistribution which would make everyone
        % better off.
    elseif strcmp(method.sim,'histogram')
        % Calculate the consumption equivalent
        % If value > 1, agents prefer policy change ( = steady state two), 
        % if 1 > value > 0, agents prefer steady state one.
        % Consumption equivalent tested against model with policy change:
        c.equivalent = ((U.two.*(1-par.sigma).*(1-par.beta)+1)...
            ./(U.one.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));
    end
end

