clear
close all
name = '40_solutions_riegler';
load(name, 'store', 'shocks');

store=store(1:24);

sol_num = size(store,2);
% sol_num =60;
store(sol_num).U = [];
U_mean = NaN(1,sol_num);
for i = 1:sol_num
    [store(i).U, U_mean(i)] = calculate_welfare( i, store, shocks, name );
end

bench = ceil(size(store,2)/2);
% mid = 30;
store(sol_num).cons_equiv = [];
cons_mean = NaN(1,sol_num);
cons_median = NaN(1,sol_num);
for i = 1:sol_num
    [store(i).cons_equiv, cons_mean(i), cons_median(i)]  = calculate_cons_equiv( i, bench, store, name );
end

em_better = NaN(1,sol_num);
unem_better = NaN(1,sol_num);
total_better = NaN(1,sol_num);
for i = 2:sol_num
    [em_better(i), unem_better(i), total_better(i)]  = per_better_two_states( i,1,store, shocks.sim_e(end,:) );
end



store(sol_num).cash_equiv = [];
cash_agg = NaN(1,sol_num);
for i = 2:sol_num
    [store(i).cash_equiv, cash_agg(i)]  = calculate_cash_equiv( i, bench, store,shocks.sim_e(end,:), name );
end

%% plot
figure(1)
title('consumption equivalent')
plot(1:sol_num, cons_median, 1:sol_num, cons_mean);
legend('median','mean')

% figure(3)
% title('Mean utility')
% plot(1:sol_num, U_mean);

figure(5)
hold on
title('Better')
plot( 1:sol_num, em_better,1:sol_num, unem_better, 1:sol_num,  total_better);
legend('Employed','Unemployed', 'Total')

figure(4)
title('Agg cash equivalent')
plot(1:sol_num, cash_agg);






