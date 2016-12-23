%% Analysis of parameters 
clear 
close 

%% Specify the parameter of interest
% 'Borrowing Constraint' 
% 'Unemployment insurance'
% 'Chance to get employed'
% etc...

%study_parameter = 'Borrowing Constraint'

k_min = [0,0.5,0.9];
T_mat = zeros(6,3);
for i =1:3
    parameters % import parameters 
    method.sim = 'histogram' 
    par.k_min = k_min(i);
    
% calling the problem solving function 
    [k,c,K,sim,store]= aiyagari_solver(par, grid, K, k, func, method, mat); 

%% Plot distribution of agents
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
    title(['kmin = ',num2str(k_min(i))]);
 
%% Saving values of Aggregates 
% calculating aggregates 
    ag1 = K.guess;
    ag2 = func.Y(K.guess);
    ag3 = func.C(K.guess);
    dev1= log(K.guess/K.rep);
    dev2= log(func.Y(K.guess)/func.Y(K.rep));
    dev3= log(func.C(K.guess)/func.C(K.rep));   
    
% saving aggregates in matrix
T_mat(:,i) = [ag1 ag2 ag3 dev1 dev2 dev3];
    
end

% save figure as pdf
saveas(gcf,'kmin.pdf')

value = {'Capital';'Output';'Consumption'};

% Create Tables 
T1 = table(T_mat(1:3,1),T_mat(4:end,1),'RowNames',value); 
T2 = table(T_mat(1:3,2),T_mat(4:end,2),'RowNames',value);
T3 = table(T_mat(1:3,3),T_mat(4:end,3),'RowNames',value);

% % Combine Tables
% 
% T_final = join(T1,T2,T3);
