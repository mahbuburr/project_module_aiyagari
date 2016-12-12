clear
close all

parameters % load paramaters
setup % load setup

method.sim = 'simulation';
% 
% tic % Set timer
% % Solve the Aiyagariy model by fixed-point iteration.
% % Function in- and output are structures.
% [k, c, K, sim, store]= aiyagari_solver(par, func, method);
% toc
% 
% % Compare welfare against frictionless model
% [c, U] = welfare_effects_rep(par, func, sim, store, K, k, method);

%% Welfare effects of two Aiyagari models with different parameters.
set(0,'DefaultFigureVisible','off'); % Do not display figures when we run the loop
tic
for i=1:2 % Get steady state values for the model with two different parameters
    if i==1 
        % e.g. par.k_min = 0.5;
        [k.one, c.one, ~, sim.one, store.one]= aiyagari_solver(par, func, method);
    elseif i==2
        % and here par.k_min = 0.25;
        par.z = 0.5;
        setup % refresh setup again for new parameter
        [k.two, c.two, ~, sim.two, store.two]= aiyagari_solver(par, func, method); % start with s.s. values from first run to speed up convergence
    end
end
toc
%set(0,'DefaultFigureVisible','on'); % Turn graphical display on again

[c, U] = welfare_effects(sim, store, k, method);