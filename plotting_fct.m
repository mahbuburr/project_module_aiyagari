function plotting_fct(par, method, gridpar, mgrid) 
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
%T_mat = zeros(12,gridpar.no,3);
T_mat = zeros(12,gridpar.no);
parameters % import parameters 

for i = 1:gridpar.no
% loop over different values of parameter you want to analyze 
    %for ii = 1:3
        if strcmp(method.analysis,'Borrowing Constraint') 
            par.k_min = par.vals(i);
        elseif strcmp(method.analysis, 'Chance to be Employed')
            par.PI_UE = par.vals(i);
        elseif strcmp(method.analysis, 'Unemployment Benefit')
            par.mu = mgrid.mu(i);
            par.PI_UE = mgrid.pi(i); % letting transition fluctuate with benefits 
        elseif strcmp(method.analysis, 'Risk Aversion')
            par.sigma = par.vals(i); 
        end
% importing functions and parameters that adapt to new input parameter value             
        setup 
% calling the problem solving function 
        [k, c, K, sim, store, mat, grid] = aiyagari_solver(par, func, method); 

% %% Plot distribution of agents
% % create subplot for each iteration
%     figure(1)
%     subplot(1,3,i)
%     if strcmp(method.sim,'simulation')
%         temp = sim.k(ceil(par.T/2):end,:);
%         hist(temp(:),100)
%         legend('number of agents')
%     elseif strcmp(method.sim,'histogram')
%         bar(grid.dist,store.distribution,'stacked')
%         legend('unemployed','employed')
%     end
% % create title matching parameter in question
%     if strcmp(method.analysis,'Borrowing Constraint') 
%         title(['kmin = ',num2str(par.k_min)]);
%     elseif strcmp(method.analysis,'Chance to be Employed')
%         title(['PI UE = ',num2str(par.PI_UE)]);
%     elseif strcmp(method.analysis,'Unemployment Benefit')
%         title(['mu = ',num2str(par.mu)]);
%     elseif strcmp(method.analysis,'Risk Aversion')
%         title(['sigma = ',num2str(par.sigma)]);
%     end
%     
%% Calculating Aggregates 
        ag1 = K.guess;
        ag2 = func.Y(K.guess);
        ag3 = func.C(K.guess); 
        w = func.w(K.guess);
        r = func.r(K.guess);
        tax = par.tau;
        labor = par.L;
        proba1 = par.PI_UE;
        proba2 = par.PI_EU;
        dev1= log(K.guess/K.rep);
        dev2= log(func.Y(K.guess)/func.Y(K.rep));
        dev3= log(func.C(K.guess)/func.C(K.rep));     
    
% saving aggregates in matrix
        T_mat(:,i) = [ag1 ag2 ag3 w r tax labor proba1 proba2 dev1 dev2 dev3]; 
        %T_mat(:,i,ii) = [ag1 ag2 ag3 w r tax labor proba1 proba2 dev1 dev2 dev3]; 
%    end
end 

%% Save figure as pdf 
saveas(gcf,'output/graph.pdf')

% %% Create Tables for Aggregates 
% 
% % set values for table
% Variable = {'Capital';'Output';'Consumption'};
% 
% % create tables and save them
% Aggregates = T_mat(1:3,1);
% log_deviations = T_mat(4:end,1);
% tables.T1 = table(Variable, Aggregates, log_deviations); 
% writetable(tables.T1,'output/Table1.txt','Delimiter',' ')
% 
% Aggregates = T_mat(1:3,2);
% log_deviations = T_mat(4:end,2);
% tables.T2 = table(Variable, Aggregates, log_deviations);
% writetable(tables.T2,'output/Table2.txt','Delimiter',' ')
% 
% Aggregates = T_mat(1:3,3);
% log_deviations = T_mat(4:end,3);
% tables.T3 = table(Variable, Aggregates, log_deviations);
% writetable(tables.T3,'output/Table3.txt','Delimiter',' ')
% T_mat_two= zeros(9,10,2);
% %% Variation of unemployment probability 
% for i =1:gridpar.no
% % loop over different values of parameter you want to analyze 
%     for k = 1:2
%         if strcmp(method.analysis,'Borrowing Constraint') 
%             par.k_min = par.vals(i);
%         elseif strcmp(method.analysis, 'Chance to be Employed')
%             par.PI_UE = par.vals(i);
%         elseif strcmp(method.analysis, 'Unemployment Benefit')
%             par.mu = mgrid.mu(i);
%             par.PI_UE = par.vals2(k); % letting transition fluctuate with benefits 
%         elseif strcmp(method.analysis, 'Risk Aversion')
%             par.sigma = par.vals(i); 
%         end
% 
% % importing functions and parameters that adapt to new input parameter value             
%         setup 
% % calling the problem solving function 
%         [k, c, K, sim, store, mat, grid] = aiyagari_solver(par, func, method); 
% 
% 
% %% Calculating Aggregates 
%         ag1 = K.guess;
%         ag2 = func.Y(K.guess);
%         ag3 = func.C(K.guess); 
%         w = func.w(K.guess);
%         r = func.r(K.guess);
%         tax = par.tau;
%         labor = par.L;
%         proba1 = par.PI_UE;
%         proba2 = par.PI_EU;
%      
% % saving aggregates in matrix
% T_mat_two(:,i,k) = [ag1 ag2 ag3 w r tax labor proba1 proba2]; 
%     end
% end 

%% Plot the grid values

figure(2)
subplot(3,3,1)
plot(mgrid.mu,T_mat(1,:))
%legend('Aggregate Capital')
xlabel('reservation wage')
ylabel('Aggregate Capital')

subplot(3,3,2)
plot(mgrid.mu,T_mat(2,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Aggregate Output')

subplot(3,3,3)
plot(mgrid.mu, T_mat(3,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Aggregate Consumption')

subplot(3,3,4)
plot(mgrid.mu,T_mat(4,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Wage')

subplot(3,3,5)
plot(mgrid.mu, T_mat(5,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Interest Rate')

subplot(3,3,6)
plot(mgrid.mu,T_mat(6,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Tax rate')

subplot(3,3,7)
plot(mgrid.mu,T_mat(7,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('employment ratio')

subplot(3,3,8)
plot(mgrid.mu,T_mat(8,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Employment probability')

subplot(3,3,9)
plot(mgrid.mu,T_mat(9,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Layoff probability')

saveas(gcf,'output/parameters.pdf')

% comparing with representative agent model
figure(3)
subplot(3,2,1)
plot(mgrid.mu,T_mat(1,:))
%legend('Aggregate Capital')
xlabel('reservation wage')
ylabel('Aggregate Capital')

subplot(3,2,2)
plot(mgrid.mu,T_mat(10,:))
%legend('Aggregate Capital')
xlabel('reservation wage')
ylabel('log-deviation Capital')

subplot(3,2,3)
plot(mgrid.mu,T_mat(2,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Aggregate Output')

subplot(3,2,4)
plot(mgrid.mu,T_mat(11,:))
%legend('Aggregate Capital')
xlabel('reservation wage')
ylabel('log-deviation Output')

subplot(3,2,5)
plot(mgrid.mu,T_mat(3,:))
%legend('Aggregate Output')
xlabel('reservation wage')
ylabel('Aggregate Consumption')

subplot(3,2,6)
plot(mgrid.mu,T_mat(12,:))
%legend('Aggregate Capital')
xlabel('reservation wage')
ylabel('log-deviation Consumption')

saveas(gcf,'output/log-deviation.pdf')
end