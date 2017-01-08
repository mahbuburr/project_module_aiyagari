%% Parameter Analysis
clear 
close 

parameters % import parameters 

% set draws for simulation
par.ind_no = 100; % number of individuals simulated
par.T = 100; % number of periods simulated

%% Specify the parameter of interest and their values
analysis = 'Risk Aversion'; % specify parameter you want to analyze: 'Borrowing Constraint' for min savings, 'Unemployment Insurance', 'Chance to be Employed','Risk Aversion'
vals = [1,2,3];
method.sim = 'simulation';

%% Create Figures for different values 
% Initiating matrix to store aggregate values    
T_mat = zeros(6,3);

for i =1:3 
% loop over different values of parameter you want to analyze 
    if strcmp(analysis,'Borrowing Constraint') 
        par.k_min = vals(i);
    elseif strcmp(analysis, 'Chance to be Employed')
        par.PI_UE = vals(i);
    elseif strcmp(analysis, 'Unemployment Insurance')
        par.mu = vals(i);
    elseif strcmp(analysis, 'Risk Aversion')
        par.sigma = vals(i); 
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
    if strcmp(analysis,'Borrowing Constraint') 
        title(['kmin = ',num2str(par.k_min)]);
    elseif strcmp(analysis,'Chance to be Employed')
        title(['PI_UE = ',num2str(par.PI_UE)]);
    elseif strcmp(analysis,'Unemployment Insurance')
        title(['mu = ',num2str(par.mu)]);
    elseif strcmp(analysis,'Risk Aversion')
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
saveas(gcf,'kmin.pdf')

%% Create Tables for Aggregates 

% set values for table
Variable = {'Capital';'Output';'Consumption'};
values = {'Aggregate variables';'log-deviation'};

% create tables and save them
Aggregates = T_mat(1:3,1);
log_deviations = T_mat(4:end,1);
T1 = table(Variable, Aggregates, log_deviations); 
writetable(T1,'Table1.txt','Delimiter',' ')

Aggregates = T_mat(1:3,2);
log_deviations = T_mat(4:end,2);
T2 = table(Variable, Aggregates, log_deviations);
writetable(T2,'Table1.txt','Delimiter',' ')

Aggregates = T_mat(1:3,3);
log_deviations = T_mat(4:end,3);
T3 = table(Variable, Aggregates, log_deviations);
writetable(T3,'Table1.txt','Delimiter',' ')
