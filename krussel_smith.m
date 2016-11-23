clear all
close all
clc

parameters; % load parameters
parameters_krussel_smith; % load parameters for Krussel - Smith

[id_shock,ag_shock]  = generate_shocks(prob,T,ind_no,U_b); % generate shocks

%% solve for general equilibrium

K_lims = [K_rep_b,grid_k(end)]; % initial limits for bisection method

% initial guesses
k_guess = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)
K_guess = (K_lims(1)+K_lims(2))/2;

d1 = 1;
iter = 0;
tic
while d1>1e-6 && iter<50 % loop for aggregate problem
    iter = iter+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Solve household problem given prices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d2 = 1;
    while d2>1e-10 % loop for household problem
             
        % calculate capital chosen next period, depending on current and
        % future employment status
        for i=1:2 % employment this period
            for j=1:2 % employment next period
                for n=1:2 % business cycle this period
                    for m=1:2 % business cycle next period
                        % capital choice next period from policy function
                        k_next(i,:,j,n,m) = interp1(grid_k,k_guess(:,j,m),k_guess(:,i,n),'linear','extrap');
%                         k_next(i,:,j) = interp1(grid_k,k_guess(:,j),k_guess(:,i),'linear','extrap'); 
                        income = mat_income(K_guess,L(n),z(n));
                        c_next(i,:,j,n,m) = max(1e-10,(1+r(K_guess,L(n),z(n))-delta)*k_guess(:,i,n) - [k_next(i,:,j,n,m)]' + income(:,m));
                    end
                end
            end
            % consumption next period from budget constraint
%             c_next(i,:,:) = max(1e-10,(1+r(K_guess)-delta)*[k_guess(:,i),k_guess(:,i)] - squeeze(k_next(i,:,:)) + mat_income(K_guess));
        end
        
%         calculate expected marginal utility of consumption next period
%         Emuc_next(:,1) = PI(1,1)*muc(c_next(1,:,1))+PI(1,2)*muc(c_next(1,:,2));  % currently unemployed in bad state
       Emuc_next(:,1,1) = prob(1,1)*muc(c_next(1,:,1,1,1)) + prob(1,2)*muc(c_next(1,:,2,1,1)) + ... % currently unemployed in bad state
                          prob(1,3)*muc(c_next(1,:,1,1,2)) + prob(1,4)*muc(c_next(1,:,2,1,2));
                       
       Emuc_next(:,2,1) = prob(2,1)*muc(c_next(2,:,1,1,1)) + prob(2,2)*muc(c_next(2,:,2,1,1)) + ... % currently employed in bad state
                          prob(2,3)*muc(c_next(2,:,1,1,2)) + prob(2,4)*muc(c_next(2,:,2,1,2));    
                     
       Emuc_next(:,1,2) = prob(3,1)*muc(c_next(1,:,1,2,1)) + prob(3,2)*muc(c_next(1,:,2,2,1)) + ... % currently unemployed in good state
                          prob(3,3)*muc(c_next(1,:,1,2,2)) + prob(3,4)*muc(c_next(1,:,2,2,2));
                      
	   Emuc_next(:,2,2) = prob(4,1)*muc(c_next(2,:,1,2,1)) + prob(4,2)*muc(c_next(2,:,2,2,1)) + ... % currently employed in good state
                          prob(4,3)*muc(c_next(2,:,1,2,2)) + prob(4,4)*muc(c_next(2,:,2,2,2));
                      
       for i=1:2
           for n=1:2
               for m=1:2
                   c_current(:,i,n) = muc_inv(beta*(1+r(K_guess,L(m),z(m))-delta)*Emuc_next(:,i,n));
                   income = mat_income(K_guess,L(n),z(n));
                   k_new(:,i,n) = (1+r(K_guess,L(n),z(n))-delta)*grid_k' + income(:,n) - c_current(:,i,n);
               end
           end
       end
                       
%         for i = 1:2
%             for n = 1:2
%                 Emuc_next(:,i,n) = prob(i,1)*muc(c_next(i,:,1,n,1)) + ... 
%                                    prob(i,2)*muc(c_next(i,:,2,n,1)) + ... 
%                                    prob(i,3)*muc(c_next(i,:,1,n,2)) + ... 
%                                    prob(i,4)*muc(c_next(i,:,2,n,2));
%                 c_current(:,i,n) = muc_inv(beta*(1+r(K_guess,L(n),z(n))-delta)*Emuc_next(:,i,n));
%             end
%         end
%         
%         for i=1:2
%             for n=1:2
                
        
%       Emuc_next(:,2) = PI(2,1)*muc(c_next(2,:,1))+PI(2,2)*muc(c_next(2,:,2));  % currently employed
        
        % calculate implied consumption this period from Euler equation
%         c_current = muc_inv(beta*(1+r(K_guess)-delta)*Emuc_next);
        
        % calculate implied capital demand from budget constraint
%         k_new = (1+r(K_guess)-delta)*mat_k + mat_income(K_guess) - c_current;
            

        % apply borrowing constraint to get new policy function
        k_new = max(k_min*K_guess,k_new); 
        
        d2_1 = norm(abs(k_new(:,1,1)-k_guess(:,1,1))./(1+abs(k_guess(:,1,1)))); % deviation between guess and new policy function
        d2_2 = norm(abs(k_new(:,2,1)-k_guess(:,2,1))./(1+abs(k_guess(:,2,1))));
        d2_3 = norm(abs(k_new(:,1,2)-k_guess(:,1,2))./(1+abs(k_guess(:,1,2))));
        d2_4 = norm(abs(k_new(:,2,2)-k_guess(:,2,2))./(1+abs(k_guess(:,2,2))));
        d2 = max([d2_1, d2_2, d2_3, d2_4]);
        
        % update policy function
        k_guess = k_guess + 0.5*(k_new-k_guess);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Find distribution of agents
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    rng('default') % reset random number generator
    sim_k = zeros(T,ind_no); % simulated values of capital stock
    sim_e = ones(T,ind_no); % simulated employment status
    sim_shock = rand(T,ind_no); % shocks for employment transition
    
    sim_e(1,1:round(L*ind_no))=2; % initial individuals that are employed
    sim_k(1,:) = K_guess; % initial capital holdings
    
    for t=2:T
        sim_e(t,sim_e(t-1,:)==1) = 1+(sim_shock(t,sim_e(t-1,:)==1)<=PI(1,2)); % new employment status of previously unemployed
        sim_e(t,sim_e(t-1,:)==2) = 1+(sim_shock(t,sim_e(t-1,:)==2)<=PI(2,2)); % new employment status of previously employed
        sim_k(t,sim_e(t,:)==1) = interp1(grid_k,k_guess(:,1),sim_k(t-1,sim_e(t,:)==1),'linear','extrap'); % capital demand of currently unemployed
        sim_k(t,sim_e(t,:)==2) = interp1(grid_k,k_guess(:,2),sim_k(t-1,sim_e(t,:)==2),'linear','extrap'); % capital demand of currently employed
    end
    
    K_demand = mean(mean(sim_k(ceil(T/2):end,:))); % average capital holdings over second half of sample
    sim_L = mean(mean(sim_e(ceil(T/2):end,:)==2)); % average employment
    
    d1 = abs(K_demand-K_guess)./(1+K_guess); % deviation between guess for capital and demanded capital stock
    
    store.K_guess(iter)=K_guess;
    store.K_next(iter)=K_demand;

    if K_demand>K_guess % update limits of bisection interval
        K_lims(1) = K_guess;
    else
        K_lims(2) = K_guess;
    end
    
    % update guess for aggregate capital stock
    K_guess = (K_lims(1)+K_lims(2))/2;
        
    disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K_guess),', K_demand: ',num2str(K_demand)])
end
toc