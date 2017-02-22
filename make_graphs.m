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
    filename = ['transitions_only' num2str(i) '.mat'];
    %filename = ['adapting_transitions_mu_parameters' num2str(i) '.mat']
    ag1(i) = load(filename, 'ag1');
    ag2(i) = load(filename, 'ag2');
    ag3(i) = load(filename, 'ag3');
    w(i) = load(filename, 'w');
    r(i) = load(filename, 'r');
    tax(i) = load(filename, 'tax');
    labor(i) = load(filename, 'labor');
    proba1(i) = load(filename, 'proba1');
    proba2(i) = load(filename, 'proba2');
    dev1(i) = load(filename, 'dev1');
    dev2(i) = load(filename, 'dev2');
    dev3(i) = load(filename, 'dev3');
    T_mat(:,i) = [ag1(i).ag1 ag2(i).ag2 ag3(i).ag3 w(i).w r(i).r tax(i).tax labor(i).labor proba1(i).proba1 proba2(i).proba2 dev1(i).dev1 dev2(i).dev2 dev3(i).dev3]; 
end


%% Plot the grid values

figure(2)
subplot(3,3,1)
plot(mgrid.mu,T_mat(1,:))
title('aggregate capital','FontWeight','bold')

subplot(3,3,2)
plot(mgrid.mu,T_mat(2,:))
title('aggregate output','FontWeight','bold')

subplot(3,3,3)
plot(mgrid.mu, T_mat(3,:))
title('aggregate consumption','FontWeight','bold')

subplot(3,3,4)
plot(mgrid.mu,T_mat(4,:))
title('wage','FontWeight','bold')

subplot(3,3,5)
plot(mgrid.mu, T_mat(5,:))
title('interest rate','FontWeight','bold')

subplot(3,3,6)
plot(mgrid.mu,T_mat(6,:))
title('tax rate','FontWeight','bold')

subplot(3,3,7)
plot(mgrid.mu,T_mat(7,:))
title('employment rate','FontWeight','bold')
%axis equal

subplot(3,3,8)
plot(mgrid.mu,T_mat(8,:))
title('job-finding probability','FontWeight','bold')

subplot(3,3,9)
plot(mgrid.mu,T_mat(9,:))
title('layoff probability','FontWeight','bold')

%mtit('Changes in the Replacement Rate', 'FontSize',14)
saveas(gcf,'output/parameters.pdf')
print('output/parameters','-dpng')


% comparing with representative agent model
figure(3)
subplot(3,2,1)
plot(mgrid.mu,T_mat(1,:))
title('aggregate capital','FontWeight','bold')

subplot(3,2,2)
plot(mgrid.mu,T_mat(10,:))
title('log-deviation capital','FontWeight','bold')

subplot(3,2,3)
plot(mgrid.mu,T_mat(2,:))
title('aggregate output','FontWeight','bold')

subplot(3,2,4)
plot(mgrid.mu,T_mat(11,:))
title('log-deviation output','FontWeight','bold')

subplot(3,2,5)
plot(mgrid.mu,T_mat(3,:))
title('aggregate consumption','FontWeight','bold')

subplot(3,2,6)
plot(mgrid.mu,T_mat(12,:))
title('log-deviation consumption','FontWeight','bold')

%mtit('Changes in the Replacement Rate','Fontsize',14)

saveas(gcf,'output/log-deviation.pdf')
print('output/log-deviation','-dpng')

% figure(4)
% subplot(2,3,1)
% plot(mgrid.mu,T_mat2(1,:,1),'b'); 
% hold on
% plot(mgrid.mu, T_mat2(1,:,2),'r'); 
% plot(mgrid.mu, T_mat2(1,:,3),'g')
% hold off 
% title('aggregate capital','FontWeight','bold')
% 
% subplot(2,3,2)
% plot(mgrid.mu,T_mat2(2,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(2,:,2),'r')
% plot(mgrid.mu, T_mat2(2,:,3),'g')
% hold off
% title('aggregate output','FontWeight','bold') 
% 
% subplot(2,3,3)
% plot(mgrid.mu, T_mat2(3,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(3,:,2),'r')
% plot(mgrid.mu, T_mat2(3,:,3),'g')
% hold off
% title('aggregate consumption','FontWeight','bold')
% legend('sigma=1','sigma=2','sigma=3')
% 
% subplot(2,3,4)
% plot(mgrid.mu,T_mat2(4,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(4,:,2),'r')
% plot(mgrid.mu, T_mat2(4,:,3),'g')
% hold off
% title('wage','FontWeight','bold')
% 
% subplot(2,3,5)
% plot(mgrid.mu, T_mat2(5,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(5,:,2),'r')
% plot(mgrid.mu, T_mat2(5,:,3),'g')
% hold off
% title('interest rate','FontWeight','bold')
% 
% subplot(2,3,6)
% plot(mgrid.mu,T_mat2(6,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(6,:,2),'r')
% plot(mgrid.mu, T_mat2(6,:,3),'g')
% hold off
% title('tax rate','FontWeight','bold')
% 
% %mtit('Changes in the Replacement Rate','Fontsize',14)
% 
% saveas(gcf,'output/parameters_sigma.pdf')
% print('output/parameters_sigma','-dpng')
% 
% %comparing with representative agent model
% figure(5)
% subplot(3,2,1)
% hold on
% plot(mgrid.mu,T_mat2(1,:,1),'b')
% plot(mgrid.mu, T_mat2(1,:,2),'r')
% plot(mgrid.mu, T_mat2(1,:,3),'g')
% title('aggregate capital','FontWeight','bold')
% 
% subplot(3,2,2)
% plot(mgrid.mu,T_mat2(10,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(10,:,2),'r')
% plot(mgrid.mu, T_mat2(10,:,3),'g')
% title('log-deviation capital','FontWeight','bold')
% legend('sigma=1','sigma=2','sigma=3')
% 
% subplot(3,2,3)
% plot(mgrid.mu,T_mat2(2,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(2,:,2),'r')
% plot(mgrid.mu, T_mat2(2,:,3),'g')
% title('aggregate output','FontWeight','bold')
% 
% 
% subplot(3,2,4)
% plot(mgrid.mu,T_mat2(11,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(11,:,2),'r')
% plot(mgrid.mu, T_mat2(11,:,3),'g')
% title('log-deviation output','FontWeight','bold')
% 
% subplot(3,2,5)
% plot(mgrid.mu,T_mat2(3,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(3,:,2),'r')
% plot(mgrid.mu, T_mat2(3,:,3),'g')
% title('aggregate consumption','FontWeight','bold')
% 
% subplot(3,2,6)
% plot(mgrid.mu,T_mat2(12,:,1),'b')
% hold on
% plot(mgrid.mu, T_mat2(12,:,2),'r')
% plot(mgrid.mu, T_mat2(12,:,3),'g')
% title('log-deviation consumption','FontWeight','bold')
% 
% %mtit('Changes in the Replacement Rate','Fontsize',14)
% saveas(gcf,'output/deviations_sigma.pdf')
% print('output/deviations_sigma','-dpng')

% figure(6)
% plot(mgrid.mu, Mat_moments(1,:),'r')
% hold on
% plot(mgrid.mu, Mat_moments(2,:),'color',[1 .5 0])
% plot(mgrid.mu, Mat_moments(5,:),'b')
% plot(mgrid.mu, Mat_moments(6,:),'c')
% xlabel('replacement rate','FontSize',14)
% ylabel('capital holdings','FontSize',14)
% mmlegend = legend('mean_{ue}','median_{ue}','mean_{e}','median_{e}');
% set(gca,'fontsize',16)
% title('Mean and Median Capital Holdings','FontSize',16)
% saveas(gcf,'output/meanandmedian.pdf')
% print('output/meanandmedian','-deps')