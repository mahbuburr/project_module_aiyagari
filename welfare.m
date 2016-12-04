clear
close all

parameters % load paramaters

tic
[k, c, K, sim, store]= f_aiyagari(par, grid, K, k, func, method, mat);
toc

%%
U = @(c) (c.^(1-sigma))./(1-sigma);
sim_c = (1+r(K_guess)-delta)*sim_k(1:end-1,:) ...
    - sim_k(2:end,:) + w(K_guess)*(1-tau)*(sim_e(1:end-1,:)-1)...
    + mu*w(K_guess)*(2-sim_e(1:end-1,:));
%life_U = sum(U(sim_c(ceil(T/2):end,:)));
figure(5)
plot(sum(sim_c,2)/ind_no)
ylabel('average consumption')

sim_u = U(sim_c);
sim_u(sim_c<0) = -Inf;