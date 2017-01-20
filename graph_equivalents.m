%% Create graph of cash equivalent for different values of unemployment benefit (mu) with mu = 0.5 as baseline
clear
close all

mu_min = 0.01;
mu_max = 0.6;
mu_n = 20;
mu = linspace(mu_min, mu_max, mu_n);
mu(7) = mu(6); % the original point does not converge
mu(6) = 0.15;
mu(8) = 0.18;

for i=1:mu_n
    filename = ['baseline_mu_' num2str(i) '.mat'];
    c(i) = load(filename, 'c');
    k(i) = load(filename, 'k');
    c_equiv(i,:) = [c(i).c.equivalent_mean, c(i).c.equivalent_median, c(i).c.equivalent_unemployed_mean, c(i).c.equivalent_unemployed_median, c(i).c.equivalent_employed_mean, c(i).c.equivalent_employed_median];
    k_equiv(i,:) = [k(i).k.equivalent_mean, k(i).k.equivalent_median, k(i).k.equivalent_unemployed_mean, k(i).k.equivalent_unemployed_median, k(i).k.equivalent_employed_mean, k(i).k.equivalent_employed_median];
end
output_baseline = 3.3539;
      
figure (1)
plot(mu', k_equiv(:,1)./ output_baseline,'r', mu', k_equiv(:,2)./ output_baseline,'--r' ...
    ,mu', k_equiv(:,3)./ output_baseline,'m', mu', k_equiv(:,5)./ output_baseline,'g' ...
    ,mu', k_equiv(:,4)./ output_baseline,'--m', mu', k_equiv(:,6)./ output_baseline,'--g')
legend('mean','median','unemployed','employed')
xlabel('unemployment benefit')
ylabel('cash equivalent / output')
refline (0,0)
 
figure (2)
plot(mu', c_equiv(:,1),'r', mu', c_equiv(:,2),'--r' ...
    ,mu', c_equiv(:,3),'m', mu', c_equiv(:,5),'g' ...
    ,mu', c_equiv(:,4),'--m', mu', c_equiv(:,6),'--g')
legend('mean','median','unemployed','employed')
xlabel('unemployment benefit')
ylabel('consumption equivalent')
refline (0,1)
 
% figure (3)
% plot(mu(2:16)', graph_var(2:16,1)./ output_baseline,'r', mu(2:16)', graph_var(2:16,2)./ output_baseline,'--r', mu(2:16)', (graph_var(2:16,3)-1).*10,'g', mu(2:16)', (graph_var(2:16,4)-1).*10,'--g')
% legend('cash equiv. mean','cash equiv. median', 'cons. equiv. mean', 'cons. equiv. median')
% xlabel('unemployment benefit')
% %ylabel('cash equivalent / output', 'consumption equivalent')
% %refline (0,0)

% %mu030= load('baseline_mu_0-3');
% mu034= load('baseline_mu_0-34');
% mu039 = load('baseline_mu_0-39');
% mu040 = load('baseline_mu_0-4');
% mu041 = load('baseline_mu_0-41');
% mu042 = load('baseline_mu_0-42');
% mu043 = load('baseline_mu_0-43');
% mu044 = load('baseline_mu_0-44');
% mu045 = load('baseline_mu_0-45');
% mu047 = load('baseline_mu_0-47');
% mu055 = load('baseline_mu_0-55');
% mu057 = load('baseline_mu_0-57');
% mu058 = load('baseline_mu_0-58');
% mu059 = load('baseline_mu_0-59');
% mu060 = load('baseline_mu_0-6');
% mu061 = load('baseline_mu_0-61');
% mu065 = load('baseline_mu_0-65');
% mu070 = load('baseline_mu_0-7');
% mu072 = load('baseline_mu_0-72');
% mu075 = load('baseline_mu_0-75');
% mu080 = load('baseline_mu_0-8');
% 
% mu = [0.34 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.47 0.5 0.55 0.57 0.58 0.59 0.6 0.61 0.65 0.7 0.72 0.75 0.8];
% graph_var = [
%             mu034.keep.k.equivalent_mean, mu034.keep.k.equivalent_median, mu034.keep.c.equivalent_mean, mu034.keep.c.equivalent_median, mu034.keep.K.two.guess;
%             mu039.keep.k.equivalent_mean, mu039.keep.k.equivalent_median, mu039.keep.c.equivalent_mean, mu039.keep.c.equivalent_median, mu039.keep.K.two.guess;
%             mu040.keep.k.equivalent_mean, mu040.keep.k.equivalent_median, mu040.keep.c.equivalent_mean, mu040.keep.c.equivalent_median, mu040.keep.K.two.guess;
%             mu041.keep.k.equivalent_mean, mu041.keep.k.equivalent_median, mu041.keep.c.equivalent_mean, mu041.keep.c.equivalent_median, mu041.keep.K.two.guess;
%             mu042.keep.k.equivalent_mean, mu042.keep.k.equivalent_median, mu042.keep.c.equivalent_mean, mu042.keep.c.equivalent_median, mu042.keep.K.two.guess;
%             mu043.keep.k.equivalent_mean, mu043.keep.k.equivalent_median, mu043.keep.c.equivalent_mean, mu043.keep.c.equivalent_median, mu043.keep.K.two.guess;
%             mu044.keep.k.equivalent_mean, mu044.keep.k.equivalent_median, mu044.keep.c.equivalent_mean, mu044.keep.c.equivalent_median, mu044.keep.K.two.guess;
%             mu045.keep.k.equivalent_mean, mu045.keep.k.equivalent_median, mu045.keep.c.equivalent_mean, mu045.keep.c.equivalent_median, mu045.keep.K.two.guess;
%             mu047.keep.k.equivalent_mean, mu047.keep.k.equivalent_median, mu047.keep.c.equivalent_mean, mu047.keep.c.equivalent_median, mu047.keep.K.two.guess;
%             0, 0, 1, 1, 25.922;
%             mu055.keep.k.equivalent_mean, mu055.keep.k.equivalent_median, mu055.keep.c.equivalent_mean, mu055.keep.c.equivalent_median, mu055.keep.K.two.guess;
%             mu057.keep.k.equivalent_mean, mu057.keep.k.equivalent_median, mu057.keep.c.equivalent_mean, mu057.keep.c.equivalent_median, mu057.keep.K.two.guess;
%             mu058.keep.k.equivalent_mean, mu058.keep.k.equivalent_median, mu058.keep.c.equivalent_mean, mu058.keep.c.equivalent_median, mu058.keep.K.two.guess;
%             mu059.keep.k.equivalent_mean, mu059.keep.k.equivalent_median, mu059.keep.c.equivalent_mean, mu059.keep.c.equivalent_median, mu059.keep.K.two.guess;
%             mu060.keep.k.equivalent_mean, mu060.keep.k.equivalent_median, mu060.keep.c.equivalent_mean, mu060.keep.c.equivalent_median, mu060.keep.K.two.guess;
%             mu061.keep.k.equivalent_mean, mu061.keep.k.equivalent_median, mu061.keep.c.equivalent_mean, mu061.keep.c.equivalent_median, mu061.keep.K.two.guess;
%             mu065.keep.k.equivalent_mean, mu065.keep.k.equivalent_median, mu065.keep.c.equivalent_mean, mu065.keep.c.equivalent_median, mu065.keep.K.two.guess;
%             mu070.keep.k.equivalent_mean, mu070.keep.k.equivalent_median, mu070.keep.c.equivalent_mean, mu070.keep.c.equivalent_median, mu070.keep.K.two.guess;
%             mu072.keep.k.equivalent_mean, mu072.keep.k.equivalent_median, mu072.keep.c.equivalent_mean, mu072.keep.c.equivalent_median, mu072.keep.K.two.guess;
%             mu075.keep.k.equivalent_mean, mu075.keep.k.equivalent_median, mu075.keep.c.equivalent_mean, mu075.keep.c.equivalent_median, mu075.keep.K.two.guess;
%             mu080.keep.k.equivalent_mean, mu080.keep.k.equivalent_median, mu080.keep.c.equivalent_mean, mu080.keep.c.equivalent_median, mu080.keep.K.two.guess];
 
%output_baseline = 2.7281; 