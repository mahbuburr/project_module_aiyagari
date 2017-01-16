function [ U ] = expec_lifetime_utility_trans( k_guess, K_demand, mu)

fixed_parameters;
tau = mu*(1-L)/L; % tax rate
mat_income = @(K) w(K)*repmat([mu,1-tau],grid_k_no,1); % matrix with income of each agent
%% ASK ABOUT CAPITAL DEMAND IN TRANSITION
c_mat = (1+r(K_demand)-delta)*mat_k-k_guess+mat_income(K_demand); 
u_mat = (c_mat.^(1-sigma)-1)/(1-sigma);
u_mat(c_mat<0) = 1e-5;

U = zeros(2,grid_k_no); %expected life time utility
dist=100;
while dist>1e-8
    EU=PI*U;
    Unew = u_mat'+betta*EU;
    dist=max(max(abs(Unew-U)));
    U=Unew;
end

end




