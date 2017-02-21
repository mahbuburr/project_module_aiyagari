%% Create graph of welfare measures for different values of unemployment benefit (mu) with mu = 0.15 as baseline
clear
close all

% Create same grid of unemployment benefits as in "welfare_analysis"
mu_min = 0.01;
mu_max = 0.6;
mu_n = 20;
mu = linspace(mu_min, mu_max, mu_n);
mu(7) = mu(6); % the original point does not converge
mu(6) = 0.15; % the original point does not converge
mu(8) = 0.18; % the original point does not converge



% Load the data
for i=1:mu_n
    %filename = ['baseline_mu_' num2str(i) '.mat'];
    filename = ['adapting_transitions_mu_' num2str(i) '.mat']
    c(i) = load(filename, 'c');
    k(i) = load(filename, 'k');
    c_equiv(i,:) = [c(i).c.equivalent_mean, c(i).c.equivalent_median, c(i).c.equivalent_unemployed_mean, c(i).c.equivalent_unemployed_median, c(i).c.equivalent_employed_mean, c(i).c.equivalent_employed_median];
    k_equiv(i,:) = [k(i).k.equivalent_mean, k(i).k.equivalent_median, k(i).k.equivalent_unemployed_mean, k(i).k.equivalent_unemployed_median, k(i).k.equivalent_employed_mean, k(i).k.equivalent_employed_median];
end

output_baseline = 3.3539;
rel_k_equiv = k_equiv./output_baseline; % Get cash equivalent relative to output
      
% Create the figures
figure (1)
plot(mu(1:17)', k_equiv(1:17,3)./ output_baseline,'m', mu(1:17)', k_equiv(1:17,5)./ output_baseline,'g' ...
    ,mu(1:17)', k_equiv(1:17,4)./ output_baseline,'--m', mu(1:17)', k_equiv(1:17,6)./ output_baseline,'--g')
legend('unemployed','employed')
xlabel('unemployment benefit')
ylabel('cash equivalent / output')
refline (0,0)
axis tight

figure (2)
plot(mu(1:13)', k_equiv(1:13,1)./ output_baseline,'r', mu(1:13)', k_equiv(1:13,2)./ output_baseline,'--r')
legend('mean','median')
xlabel('unemployment benefit')
ylabel('cash equivalent / output')
refline (0,0)


figure (3)
plot(mu', c_equiv(:,3),'m', mu', c_equiv(:,5),'g' ...
    ,mu', c_equiv(:,4),'--m', mu', c_equiv(:,6),'--g')
legend('unemployed','employed')
xlabel('unemployment benefit')
ylabel('consumption equivalent')
refline (0,1)
axis tight

figure (4)
plot(mu', c_equiv(:,1),'r', mu', c_equiv(:,2),'--r')
legend('mean','median')
xlabel('unemployment benefit')
ylabel('consumption equivalent')
refline (0,1)
axis tight

figure (5)
plot(mu(1:9)', k_equiv(1:9,1)./ output_baseline,'r', mu(1:9)', (c_equiv(1:9,1)-1).*100, 'g' ...
    ,mu(1:9)', k_equiv(1:9,2)./ output_baseline,'--r' , mu(1:9)', (c_equiv(1:9,2)-1).*100,'--g')
legend('cash equivalent', 'consumption equivivalent')
xlabel('unemployment benefit')
%ylabel('cash equivalent / output', 'consumption equivalent')
refline (0,0)
axis tight