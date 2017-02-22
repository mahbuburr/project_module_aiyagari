% plot the graphs
clear 
close all 

% Create grid of unemployment benefits for the analysis
gridpar.start = 0.01;
gridpar.end = 0.6;
gridpar.no = 20;
mgrid.mu = linspace(gridpar.start, gridpar.end, gridpar.no);
mgrid.mu(7) = mgrid.mu(6); % the original point does not converge
mgrid.mu(6) = 0.15; % the original point does not converge
mgrid.mu(8) = 0.18; % the original point does not converge



for i=1:gridpar.no
    %filename = ['baseline_mu_' num2str(i) '.mat'];
    filename = ['adapting_transitions_mu_parameters' num2str(i) '.mat']
    T_mat(:,i) = [ag1 ag2 ag3 w r tax labor proba1 proba2 dev1 dev2 dev3];
end



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