function [Mat_moments] = plotting_fct(par, method, gridpar, mgrid) 
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
T_mat = zeros(12,gridpar.no);
Mat_moments = zeros(8,gridpar.no);

T_mat2 = zeros(12,gridpar.no,3);
Mat_moments2 = zeros(2,4,3,3);

parameters % import parameters 

for i = 11%:gridpar.no
tic
iteration = i
% loop over different values of parameter you want to analyze 
     %for ii = 1:3
     
        if strcmp(method.analysis,'Borrowing Constraint') 
            par.k_min = par.vals(i);
        elseif strcmp(method.analysis, 'Chance to be Employed')
            par.PI_UE = mgrid.pi(i);
        elseif strcmp(method.analysis, 'Unemployment Benefit')
            par.mu = mgrid.mu(i);
%               par.sigma = par.vals2(ii);
           % par.PI_UE = mgrid.pi(i); % letting transition fluctuate with benefits 
        elseif strcmp(method.analysis, 'Risk Aversion')
            par.sigma = par.vals(i); 
        end
% importing functions and parameters that adapt to new input parameter value             
        setup 
%method.HH = 'FP'; % Depending on the mu, you might have to change it to 'FP' or 'FPend' to converge
%method.agg = 'bisection'; % Depending on the mu, you might have to change it to 'bisection' or 'bisectio' to converge
% calling the problem solving function 
        [k, c, K, sim, store, mat, grid] = aiyagari_solver(par, func, method); 

% %% Plot distribution of agents
    figure(1)
    if strcmp(method.sim,'simulation')
        temp = sim.k(ceil(par.T/2):end,:);
        hist(temp(:),100)
        legend('number of agents')
    elseif strcmp(method.sim,'histogram')
        bar(grid.dist,store.distribution,'stacked')
        legend('unemployed','employed')
        %set(gca,'fontsize',25)
        
    end

% print(['output/distribution=  ', num2str(i)],'-dpng')
% 
% 
% end
%% Calculate and save the moments 


capital_grid = [transpose(grid.dist) transpose(grid.dist)];
distr = store.distribution.*capital_grid;

% calculating the mean [unemployed employed]
% note that you need to rescale the probabilities!!!
mean_dis = sum(distr).*[1/(1-par.L) 1/par.L];

sum1 = cumsum(store.distribution(:,1));
sum2 = cumsum(store.distribution(:,2));
[~, index1]=min(abs(sum1-((sum1(end))/2)));
[~, index2]=min(abs(sum2-((sum2(end))/2)));
closestValues=sum1(index1);% to check value for unemployed
closestValues2 =sum2(index2); % to check value for employed

kept.mean_unemployed = mean_dis(1);
kept.mean_employed = mean_dis(2);
kept.median_unemployed = capital_grid(index1); %using index to get median for unemployed
kept.median_employed = capital_grid(index2);% using index to get median for employed
kept.std_unemployed =  sqrt(mean((grid.dist - mean_dis(1)).^2)); %the standard deviation of the unemployed
kept.std_employed =  sqrt(mean((grid.dist - mean_dis(2)).^2)); % standard deviation of the employed 
kept.skewness_unemployed = mean(((grid.dist - mean_dis(1))/kept.std_unemployed).^3); %skewness of the unemployed
kept.skewness_employed = mean(((grid.dist - mean_dis(2))/kept.std_employed).^3); % skewness of the employed 

% store moments
%Mat_moments = [mean_dis(1) median_dis(1) std_dis1 skewness_dis1 mean_dis(2) median_dis(2) std_dis2 skewness_dis2];
savename = ['moments_transitions_only' num2str(i) '.mat'];
save(savename, '-struct','kept'); % Save the values

% Mat_moments2(:,:,i,ii) = [mean_dis(1) median_dis(1) std_dis1 skewness_dis1;
%                           mean_dis(2) median_dis(2) std_dis2 skewness_dis2];
%               
% %% Calculating Aggregates
%         keep.ag1 = K.guess;
%         keep.ag2 = func.Y(K.guess);
%         keep.ag3 = func.C(K.guess); 
%         keep.w = func.w(K.guess);
%         keep.r = func.r(K.guess);
%         keep.tax = par.tau;
%         keep.labor = par.L;
%         keep.proba1 = par.PI_UE;
%         keep.proba2 = par.PI_EU;
%         keep.dev1= log(K.guess/K.rep);
%         keep.dev2= log(func.Y(K.guess)/func.Y(K.rep));
%         keep.dev3= log(func.C(K.guess)/func.C(K.rep));     
%     
% % saving the values
%         filename = ['transitions_only' num2str(i) '.mat']
%         %filename = ['adapting_transitions_mu_parameters' num2str(i) '.mat'];
%         %filename = ['baseline_mu_parameters' num2str(i) '.mat'];
%         save (filename, '-struct','keep'); % Save the values
%         %T_mat2(:,i,ii) = [ag1 ag2 ag3 w r tax labor proba1 proba2 dev1 dev2 dev3]; 
%     % end
    toc
end 

%% Save figure as pdf 
% saveas(gcf,'output/distribution.pdf')
% print('output/distribution','-dpng')


end