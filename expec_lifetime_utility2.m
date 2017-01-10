period = 2500;
e = sim_e(period,:); 
c_mat = (1+r(K_demand)-delta)*mat_k-k_guess+mat_income(K_demand); 
u_mat = (c_mat.^(1-sigma)-1)/(1-sigma);
u_mat(c_mat<0) = 1e-5;

U = zeros(2,grid_k_no); %expected life time utility
dist=100;
iter1 = 0;
while dist>1e-8
    EU=PI*U;
    Unew = u_mat'+beta*EU;
    dist=max(max(abs(Unew-U)));
    U=Unew;
end


store.U = NaN(1,size(e,2));
for individual = 1:size(e,2)
    store.U(individual) = interp1(grid_k,U(e(individual),:),sim_k(period,individual),'linear','extrap');
end




