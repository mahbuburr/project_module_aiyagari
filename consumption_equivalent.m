function [ c ] = consumption_equivalent (par, method, U)
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
        c.equivakent_median = median(median(c.equivalent));
        
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

