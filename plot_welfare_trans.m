clear
close all
name = '20_solutions_riegler';
load(name, 'store', 'shocks');

sol_num = size(store,2);
store(sol_num).U = [];
U_mean = NaN(1,sol_num);
for i = 1:sol_num
    [store(i).U, U_mean(i)] = calculate_welfare( i, store, shocks, name );
end

% mid = ceil(size(store,2)/2);
mid = 2;
store(sol_num).cons_equiv = [];
cons_mean = NaN(1,sol_num);
cons_median = NaN(1,sol_num);
for i = 1:sol_num
    [store(i).cons_equiv, cons_mean(i), cons_median(i)]  = calculate_cons_equiv( i, mid, store, name );
end

store(sol_num).cash_equiv = [];
cash_agg = NaN(1,sol_num);
for i = 1:sol_num
    [store(i).cash_equiv, cash_agg(i)]  = calculate_cons_equiv( i, mid, store, name );
end

%% plot
figure(1)
plot(1:sol_num, cons_median);

figure(2)
plot(1:sol_num, cons_mean);

figure(3)
plot(1:sol_num, U_mean);

figure(4)
plot(1:sol_num, cash_agg);





