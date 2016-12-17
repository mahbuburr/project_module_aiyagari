clear all
close all

warning('off','all')
parameters; % load parameters
[prob, probaux] = probability_matrix(U_b, PI_UE_b, U_g, PI_UE_g, B_p, PI_bg, grid_k_no, grid_K_no, ag_states_no, id_states_no);
[id_shock,ag_shock]  = generate_shocks(prob,T,ind_no,U_b); % generate shocks
% load('shocks_vector_internet.mat');
% id_shock = idshock;
% ag_shock = agshock;



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

sim_k = zeros(T,ind_no); % simulated values of capital stock
sim_k(1,:) = kss; % initial capital holdings

dif_B = 1;
iter = 0;
tic
while dif_B>1e-8 %&& iter<50 % loop for aggregate problem
    iter = iter+1
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
        
        a = load('individual_from_int.mat');
        % Bad aggregate state and unemployed idiosyncratic state
        % future capital state k''
        k_next_bu=interpn(grid_k,grid_K,k_guess(:,:,1,1),k_guess,K_guess,'cubic');
        % future consumption (c')
        c_next_bu=max(1e-10,r(K_guess,L(1),z(1)).*k_guess+mat_income(K_guess,L(1),z(1),e(1))+(1-delta).*k_guess-k_next_bu);
        % marginal utility of future consumption
        mu_next_bu=muc_inv(c_next_bu); 
        
        % Bad aggregate state and employed idiosyncratic state
        k_next_be=interpn(grid_k,grid_K,k_guess(:,:,1,2),k_guess,K_guess,'cubic');
        c_next_be=max(1e-10,r(K_guess,L(1),z(1)).*k_guess+mat_income(K_guess,L(1),z(1),e(2))+(1-delta).*k_guess-k_next_be);
        mu_next_be=muc_inv(c_next_be);
        
        % Good aggregate state and unemployed idiosyncratic state
        k_next_gu=interpn(grid_k,grid_K,k_guess(:,:,2,1),k_guess,K_guess,'cubic');
        c_next_gu=max(1e-10,r(K_guess,L(2),z(2)).*k_guess+mat_income(K_guess,L(2),z(2),e(1))+(1-delta).*k_guess-k_next_gu);
        mu_next_gu=muc_inv(c_next_gu);
        
        % Good aggregate state and employed idiosyncratic state
        k_next_ge=interpn(grid_k,grid_K,k_guess(:,:,2,2),k_guess,K_guess,'cubic');
        c_next_ge=max(1e-10,r(K_guess,L(2),z(2)).*k_guess+mat_income(K_guess,L(2),z(2),e(2))+(1-delta).*k_guess-k_next_ge);
        mu_next_ge=muc_inv(c_next_ge);
        
        
        % calculate expected marginal utility of consumption next period
        Emuc_next =(mu_next_bu.*(1-delta+r(K_guess,L(1),z(1)))).*probaux(:,:,:,:,1) + ...
            (mu_next_be.*(1-delta+r(K_guess,L(1),z(1)))).*probaux(:,:,:,:,2) + ...
            (mu_next_gu.*(1-delta+r(K_guess,L(2),z(2)))).*probaux(:,:,:,:,3) + ...
            (mu_next_ge.*(1-delta+r(K_guess,L(2),z(2)))).*probaux(:,:,:,:,4);
        
        c_current = muc_inv(beta*Emuc_next);
        k_new = NaN(100,4,2,2);
        k_new(:,:,1,1) = wealth(grid_K, grid_k, L(1), z(1), e(1)) - c_current(:,:,1,1);
        k_new(:,:,1,2) = wealth(grid_K, grid_k, L(1), z(1), e(2)) - c_current(:,:,1,2);
        k_new(:,:,2,1) = wealth(grid_K, grid_k, L(2), z(2), e(1)) - c_current(:,:,2,1);
        k_new(:,:,2,2) = wealth(grid_K, grid_k, L(2), z(2), e(2)) - c_current(:,:,2,2);
        k_new=min(max(k_min, k_new),k_max); % apply borrowing constraint to get new policy function

        d2=max(max(max(max(abs(k_new-k_guess)))));
    
        % update policy function
        k_guess = k_guess + 0.7*(k_new-k_guess);
    end
%     load('policy.mat', 'k_guess');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Find distribution of agents
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Solving aggregate problem');

    K_demand = zeros(T,1);
    for t = 1:T
        K_demand(t) = mean(sim_k(t,:));   
        K_demand(t) = min(max(K_min, K_demand(t)),K_max);
        
        k_aux = squeeze(interpn(grid_k, grid_K, z_s, e_s, k_guess, grid_k, K_demand(t), ag_shock(t), e_s,'cubic'));
        sim_k(t+1,:) = interpn(grid_k,e_s,k_aux,sim_k(t,:),id_shock(t,:),'cubic'); 
        sim_k(t+1,:) = min(max(k_min, sim_k(t+1,:)),k_max);
    end
    b=load('agg_from_int.mat');
    disp('Regression');
    ibad=0;           % count how many times the aggregate shock was bad
    igood=0;          % count how many times the aggregate shock was good
    xbad=0;  ybad=0;  % regression-variables for a bad state
    xgood=0; ygood=0; % regression-variables for a good state
    for i=100+1:T-1
        if ag_shock(i)==1
            ibad=ibad+1;
            xbad(ibad,1)=log(K_demand(i));
            ybad(ibad,1)=log(K_demand(i+1));
        else
            igood=igood+1;
            xgood(igood,1)=log(K_demand(i));
            ygood(igood,1)=log(K_demand(i+1));
        end
    end
    
    [B1(1:2),s2,s3,s4,s5]=regress(ybad,[ones(ibad,1) xbad]);R2bad=s5(1);
    % run the OLS regression ln(km')=B(1)+B(2)*ln(km) for a bad agg. state
    % and compute R^2 (which is the first statistic in s5)
    [B1(3:4),s2,s3,s4,s5]=regress(ygood,[ones(igood,1) xgood]);R2good=s5(1);
    % make the OLS regression ln(km')=B(3)+B(4)*ln(km) for a good agg. state
    % and compute R^2 (which is the first statistic in s5)
    
    dif_B=norm(B-B1) % compute the difference between the initial and obtained
    % vector of coefficients
    
    % To ensure that initial capital distribution comes from the ergodic set,
    % we use the terminal distribution of the current iteration as initial
    % distribution for a subsequent iteration. When the solution is sufficiently
    % accurate, dif_B<(criter_B*100), we stop such an updating and hold the
    % distribution "kcross" fixed for the rest of iterations. ·
    
    if dif_B>(1e-8*100)
        sim_k(1,:)=sim_k(T,:); % the new capital distribution  replaces the old one
    end
    
    B=B1*update_B+B*(1-update_B); % update the vector of the ALM coefficients
    % according to the rule (9) in the paper

    
end
toc
j = load('Solution_to_model');