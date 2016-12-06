function [ output_args ] = welfare_effects( input_args )
%WELFARE EFFECTS Calculates whether agents prefer a policy change or not
%   Detailed explanation goes here

% Define function to calculate utility of consumption.
func.U = @(c) (c.^(1-par.sigma))./(1-par.sigma);

if strcmp(method.sim ,'simulation')

    ind_no = size(sim.k,2);
    % Return consumtion values for all agents.
    sim.c = (1+func.r(K.guess)-par.delta)*sim.k(1:end-1,:) ...
        - sim.k(2:end,:) + func.w(K.guess)*(1-par.tau)*(sim.e(1:end-1,:)-1)...
        + par.mu*func.w(K.guess)*(2-sim.e(1:end-1,:));

    sim.u = func.U(sim.c);
    sim.u(sim.c<0) = -Inf; % Make sure consumtion is always positive.
    % Calculate the life time utility when simulation has converged (about half
    % the size of the simulation).
    lifetime_U = par.beta.^(1:ceil((ind_no-1)/2))*sim.u(ceil((ind_no-1)/2):end,:);
    lifetime_U_rep = sum(par.beta.^(1:floor((ind_no-1)/2))*func.U(func.C(K.guess)));
    
elseif strcmp(method.sim , 'histogram')
% still need to figure this out


%% Calculate the consumption equivalent
% If value>1, agents prefer policy change, if 1>value>0, agents prefer
% old model.
% Consumption equivalent tested against frictionless benchmark model.
c.equivalent_bench = ((lifetime_U.*(1-par.sigma).*(1-par.beta)+1)...
    ./(lifetime_U_rep*(1-par.sigma)*(1-par.beta)+1)).^(1/(1-par.sigma)); 
histogram(c.equivalent_bench)
% Calculate average and median of consumption equivalent. If this aggregate
% >1, there exists a (lump-sum) redistribution which would make everyone
% better off.
c.equivalent_bench_mean = mean(c.equivalent_bench);
c.equivalent_bench_median = median(c.equivalent_bench);

end

