period = 2500;
e = sim_e(period,:); %employment status 0 or 1
e(e==1)=0; 
e(e==2)=1;
c_initial = (1+r(K_demand)-delta)*sim_k(period,:)-sim_k(period+1,:)+w(K_demand)*(1-tau)*(e) + mu*w(K_demand)*(1-e); %consumption for all individuals 
u_initial = (c_initial.^(1-sigma)-1)/(1-sigma); %utility from cons. 

c_inital = 

store.U = NaN(1,size(e,2));
for individual = 1:size(e,2)
    U = zeros(2,grid_k_no); %expected life time utility
    dist=100;
    while dist>1e-8
        EU=PI*U;
        Unew = u_initial(individual)+beta*EU;
        dist=max(max(abs(Unew-U)));
        U=Unew;
    end
    store.U(individual) = interp1(grid_k,U(e(individual)+1,:),sim_k(period,individual),'linear','extrap');
end





