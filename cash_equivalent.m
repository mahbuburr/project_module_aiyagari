function [ k ] = cash_equivalent( par, method, grid, sim, U )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(method.sim ,'simulation')
        ind_no = size(sim.one.k,2);
        T = size(sim.one.k,1);
        k.equivalent = NaN(ceil(T/2),ind_no);
        k.compensated = NaN(ceil(T/2),ind_no);
        tic
        for t = ceil((T+1)/2):T
            k.compensated(t-2500,sim.one.e(t,:)==1) = interp1(U.two.lifetime(1,:), grid.one.k, U.one.one(t-2500,sim.one.e(t,:)==1), 'linear', 'extrap');
            k.compensated(t-2500,sim.one.e(t,:)==2) = interp1(U.two.lifetime(2,:), grid.one.k, U.one.one(t-2500,sim.one.e(t,:)==2), 'linear', 'extrap');
            k.equivalent(t-2500,:) = k.compensated(t-2500,:) - sim.one.k(t,:);   
        end
        toc
    elseif strcmp(method.sim, 'histogram')
        % to be done
end

