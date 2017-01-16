clear
close all

load('10solutions_riegler.mat', 'store', 'shocks');
store(10).U = [];
U_mean = NaN(1,10);
for i = 1:10
    [store(i).U, U_mean(i)] = calculate_welfare( i, store, shocks );
end

mid = ceil(size(store,2)/2);
store(10).cons_equiv = [];
cons_mean = NaN(1,10);
cons_median = NaN(1,10);
for i = 1:10
    [store(i).cons_equiv, cons_mean(i), cons_median(i)]  = calculate_cons_equiv( i, mid, store );
end

store(10).cash_equiv = [];
cash_agg = NaN(1,10);
for i = 1:10
    [store(i).cash_equiv, cash_agg(i)]  = calculate_cons_equiv( i, mid, store );
end

%% plot
figure(1)
plot(1:10, cons_median);

figure(2)
plot(1:10, cons_mean);

figure(3)
plot(1:10, U_mean);

figure(4)
plot(1:10, cash_agg);





