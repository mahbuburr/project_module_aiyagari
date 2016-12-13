%% Analysis of parameters 
clear 
close 

parameters % import parameters 
method.sim = 'histogram' 
%% Specify the parameter of interest
% 'Borrowing Constraint' 
% 'Unemployment insurance'
% 'Chance to get employed'
% etc...

%study_parameter = 'Borrowing Constraint'

k_min = [0,0.5,0.9];
for i =1:3
    par.k_min = k_min(i);
% calling the problem solving function 
    [k,c,K,sim,store]= aiyagari_solver(par, grid, K, k, func, method, mat); 

% Plot distribution of agents

    figure(1)
    
    subplot(1,3,i)
    if strcmp(method.sim,'simulation')
        temp = sim.k(ceil(T/2):end,:);
        hist(temp(:),100)
        legend('number of agents')
    elseif strcmp(method.sim,'histogram')
        bar(grid.dist,store.distribution,'stacked')
        legend('unemployed','employed')
    end

end