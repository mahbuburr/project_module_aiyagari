clear all
close all

parameters; % load parameters
[id_shock,ag_shock]  = generate_shocks(prob,T,ind_no,U_b); % generate shocks

%% solve for general equilibrium
% initial guesses

% next-period individual capital (k') depends on four state variables: 
% individual k, aggregate k, aggregate shock, idiosyncratic shock 
k_guess=zeros(grid_k_no,grid_K_no,ag_states_no,id_states_no); 

% Initial policy function
for i=1:4
   for j=1:2
      for h=1:2
         k_guess(:,i,j,h)=0.9*grid_k;
      end
   end
end

K_guess=zeros(grid_k_no,grid_K_no,ag_states_no,id_states_no);

dif_B = 1;
iter = 0;
tic
while dif_B>1e-6 && iter<50 % loop for aggregate problem
    iter = iter+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Solve household problem given prices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Solving individual problem');
    d2 = 1;
    while d2>1e-8 % loop for household problem
        
        K_guess(:,:,1,1)=exp(repmat(B(1),grid_k_no,grid_K_no)+B(2)*log(repmat(grid_K',100,1)));
        K_guess(:,:,1,2)=exp(repmat(B(1),grid_k_no,grid_K_no)+B(2)*log(repmat(grid_K',100,1)));
        K_guess(:,:,2,1)=exp(repmat(B(3),grid_k_no,grid_K_no)+B(4)*log(repmat(grid_K',100,1)));
        K_guess(:,:,2,2)=exp(repmat(B(3),grid_k_no,grid_K_no)+B(4)*log(repmat(grid_K',100,1)));
        K_guess=min(max(K_min, K_guess),K_max); % restricting K' to be in [K_min,K_max] range
        
        
        % Bad aggregate state and unemployed idiosyncratic state
        % future capital state k''
        k_next_bu=interpn(grid_k,grid_K,k_guess(:,:,1,1),k_guess,K_guess,'spline');
        % future consumption (c')
        c_next_bu=max(1e-10,r_mat(grid_K,L(1),z(1)).*k_guess+mat_income(grid_K,L(1),z(1),e(1))-k_next_bu);
        % marginal utility of future consumption
        mu_next_bu=muc_inv(c_next_bu); 
        
        % Bad aggregate state and employed idiosyncratic state
        k_next_be=interpn(grid_k,grid_K,k_guess(:,:,1,2),k_guess,K_guess,'spline');
        c_next_be=max(1e-10,r_mat(grid_K,L(1),z(1)).*k_guess+mat_income(grid_K,L(1),z(1),e(2)) - k_next_be);
        mu_next_be=muc_inv(c_next_be);
        
        % Good aggregate state and unemployed idiosyncratic state
        k_next_gu=interpn(grid_k,grid_K,k_guess(:,:,2,1),k_guess,K_guess,'spline');
        c_next_gu=max(1e-10,r_mat(grid_K,L(2),z(2)).*k_guess+mat_income(grid_K,L(2),z(2),e(1))-k_next_gu);
        mu_next_gu=muc_inv(c_next_gu);
        
        % Good aggregate state and employed idiosyncratic state
        k_next_ge=interpn(grid_k,grid_K,k_guess(:,:,2,2),k_guess,K_guess,'spline');
        c_next_ge=max(1e-10,r_mat(grid_K,L(2),z(2)).*k_guess+mat_income(grid_K,L(2),z(2),e(2)) - k_next_ge);
        mu_next_ge=muc_inv(c_next_ge);
        
        
        % calculate expected marginal utility of consumption next period
        Emuc_next =(mu_next_bu.*(r_mat(grid_K,L(1),z(1)))).*probaux(:,:,:,:,1) + ...
            (mu_next_be.*r_mat(grid_K,L(1),z(1))).*probaux(:,:,:,:,2) + ...
            (mu_next_gu.*r_mat(grid_K,L(1),z(1))).*probaux(:,:,:,:,3) + ...
            (mu_next_ge.*r_mat(grid_K,L(1),z(1))).*probaux(:,:,:,:,4);
        
        c_current = muc_inv(beta*Emuc_next);
        k_new = NaN(100,4,2,2);
        k_new(:,:,1,1) = (1-delta+repmat(r(grid_K',L(1),z(1)),grid_k_no,1)).*repmat(grid_k,1,grid_K_no)+repmat(income(grid_K,L(1),z(1),e(1)),grid_k_no,1)-c_current(:,:,1,1);
        k_new(:,:,1,2) = (1-delta+repmat(r(grid_K',L(1),z(1)),grid_k_no,1)).*repmat(grid_k,1,grid_K_no)+repmat(income(grid_K,L(1),z(1),e(2)),grid_k_no,1)-c_current(:,:,1,2);
        k_new(:,:,2,1) = (1-delta+repmat(r(grid_K',L(2),z(2)),grid_k_no,1)).*repmat(grid_k,1,grid_K_no)+repmat(income(grid_K,L(2),z(2),e(1)),grid_k_no,1)-c_current(:,:,2,1);
        k_new(:,:,2,2) = (1-delta+repmat(r(grid_K',L(2),z(2)),grid_k_no,1)).*repmat(grid_k,1,grid_K_no)+repmat(income(grid_K,L(2),z(2),e(2)),grid_k_no,1)-c_current(:,:,2,2);
        k_new=min(max(k_min, k_new),k_max); % apply borrowing constraint to get new policy function

        d2=max(max(max(max(abs(k_new-k_guess)))));
    
        % update policy function
        k_guess = k_guess + 0.7*(k_new-k_guess);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Find distribution of agents
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Solving aggregate problem');
    sim_k = zeros(T,ind_no); % simulated values of capital stock
    kcross=zeros(1,N)+kss;
    sim_k(1,:) = K_guess; % initial capital holdings
    
    for t=2:T
            if ag_shock(t) == 1
                 sim_k(t,id_shock(t,:)==1) = interp1(grid_k,k_guess(:,1,1),sim_k(t-1,id_shock(t,:)==1),'linear','extrap'); % capital demand of currently unemployed
                 sim_k(t,id_shock(t,:)==2) = interp1(grid_k,k_guess(:,2,1),sim_k(t-1,id_shock(t,:)==2),'linear','extrap'); % capital demand of currently employed
            else
                sim_k(t,id_shock(t,:)==1) = interp1(grid_k,k_guess(:,1,2),sim_k(t-1,id_shock(t,:)==1),'linear','extrap'); % capital demand of currently unemployed
                sim_k(t,id_shock(t,:)==2) = interp1(grid_k,k_guess(:,2,2),sim_k(t-1,id_shock(t,:)==2),'linear','extrap'); % capital demand of cu
            end
    end
 
end
toc