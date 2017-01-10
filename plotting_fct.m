function plotting_fct(par, method) 
%% PLOTTING FUNCTION: The function creates a figure with three different
% subplots for a vector of a parmater chosen for analysis. For each value
% the parameter takes, the function creates a table containing the
% aggregate values. Both, the three tables and the figure are saved in an
% output folder. 
%
% Input variables:
%   par = vector of parameters:
%       ind_no = number of individuals simulated
%       T = number of periods simulated
%       vals = vector of three parameters used to run the analysis 
%   method = method of simulation choice and method of analysis choice
%
%% Create Figures for different values 
% Initiating matrix to store aggregate values    
T_mat = zeros(6,3);

parameters % import parameters 

for i =1:3 
% loop over different values of parameter you want to analyze 
    if strcmp(method.analysis,'Borrowing Constraint') 
        par.k_min = par.vals(i);
    elseif strcmp(method.analysis, 'Chance to be Employed')
        par.PI_UE = par.vals(i);
    elseif strcmp(method.analysis, 'Unemployment Benefit')
        par.mu = par.vals(i);
    elseif strcmp(method.analysis, 'Risk Aversion')
        par.sigma = par.vals(i); 
    end
% importing functions and parameters that adapt to new input parameter value             
    setup 
% calling the problem solving function 
    [k, c, K, sim, store, mat, grid] = aiyagari_solver(par, func, method); 

%% Plot distribution of agents
% create subplot for each iteration
    figure(1)
    subplot(1,3,i)
    if strcmp(method.sim,'simulation')
        temp = sim.k(ceil(par.T/2):end,:);
        hist(temp(:),100)
        legend('number of agents')
    elseif strcmp(method.sim,'histogram')
        bar(grid.dist,store.distribution,'stacked')
        legend('unemployed','employed')
    end
% create title matching parameter in question
    if strcmp(method.analysis,'Borrowing Constraint') 
        title(['kmin = ',num2str(par.k_min)]);
    elseif strcmp(method.analysis,'Chance to be Employed')
        title(['PI UE = ',num2str(par.PI_UE)]);
    elseif strcmp(method.analysis,'Unemployment Benefit')
        title(['mu = ',num2str(par.mu)]);
    elseif strcmp(method.analysis,'Risk Aversion')
        title(['sigma = ',num2str(par.sigma)]);
    end
    
%% Calculating Aggregates 
    ag1 = K.guess;
    ag2 = func.Y(K.guess);
    ag3 = func.C(K.guess);
    dev1= log(K.guess/K.rep);
    dev2= log(func.Y(K.guess)/func.Y(K.rep));
    dev3= log(func.C(K.guess)/func.C(K.rep));   
    
% saving aggregates in matrix
T_mat(:,i) = [ag1 ag2 ag3 dev1 dev2 dev3]; 
end 

%% Save figure as pdf 
saveas(gcf,'output/graph.pdf')

%% Create Tables for Aggregates 

% set values for table
Variable = {'Capital';'Output';'Consumption'};

% create tables and save them
Aggregates = T_mat(1:3,1);
log_deviations = T_mat(4:end,1);
tables.T1 = table(Variable, Aggregates, log_deviations); 
writetable(tables.T1,'output/Table1.txt','Delimiter',' ')

Aggregates = T_mat(1:3,2);
log_deviations = T_mat(4:end,2);
tables.T2 = table(Variable, Aggregates, log_deviations);
writetable(tables.T2,'output/Table2.txt','Delimiter',' ')

Aggregates = T_mat(1:3,3);
log_deviations = T_mat(4:end,3);
tables.T3 = table(Variable, Aggregates, log_deviations);
writetable(tables.T3,'output/Table3.txt','Delimiter',' ')
end
