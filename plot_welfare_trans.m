clear
close all

load('10solutions_riegler_trans.mat', 'store', 'shocks');
store(11).U = [];
for n = 2:11
    store(n).U = calculate_welfare_trans( n, store, shocks );
end


