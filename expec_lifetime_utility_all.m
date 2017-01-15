function [ store ] = expec_lifetime_utility_all( k_guess, K_demand, mu, sim_e, sim_k )

fixed_parameters;
tau = mu*(1-L)/L; % tax rate
mat_income = @(K) w(K)*repmat([mu,1-tau],grid_k_no,1); % matrix with income of each agent

period = 2500;
e = sim_e(period,:); 
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


store.U = NaN(1,size(e,2));
for individual = 1:size(e,2)
    Uaux = NaN(1,2500);
    for t = period:T
        Uaux(t) = interp1(grid_k,U(e(individual),:),sim_k(period,individual),'linear','extrap');
    end
    store.U(individual) = mean(Uaux);
end
end



