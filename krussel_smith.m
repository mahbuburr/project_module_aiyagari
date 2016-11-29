clear all
close all

parameters; % load parameters
% parameters_krussel_smith; % load parameters for Krussel - Smith

[id_shock,ag_shock]  = generate_shocks(prob,T,ind_no,U_b); % generate shocks

%% solve for general equilibrium

K_lims = [K_rep_b,grid_k(end)]; % initial limits for bisection method

% initial guesses
% k_guess = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)
% k_guess = repmat(0.9*37.98,[100,4,2,2]);
K_guess = (K_lims(1)+K_lims(2))/2;

k_guess=zeros(100,4,2,2); % next-period individual 
   % capital (k') depends on four state variables: individual k, aggregate k, 
   % aggregate shock, idiosyncratic shock 

% Initial capital function

for i=1:4
   for j=1:2
      for h=1:2
         k_guess(:,i,j,h)=0.9*grid_k;
      end
   end
end

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
        kmprime=zeros(grid_k_no,ngridkm,2,2);
        kmprime(:,:,1,1)=exp(B(1)*ones(grid_k_no,ngridkm)+B(2)*log(kmaux(:,:,1,1)));
        kmprime(:,:,1,2)=exp(B(1)*ones(grid_k_no,ngridkm)+B(2)*log(kmaux(:,:,1,2)));
        kmprime(:,:,2,1)=exp(B(3)*ones(grid_k_no,ngridkm)+B(4)*log(kmaux(:,:,2,1)));
        kmprime(:,:,2,2)=exp(B(3)*ones(grid_k_no,ngridkm)+B(4)*log(kmaux(:,:,2,2)));
        kmprime=(kmprime>=km_min).*(kmprime<=km_max).*kmprime+(kmprime<km_min)*km_min+(kmprime>km_max)*km_max;
        % calculate capital chosen next period, depending on current and
        % future employment status
%         for i=1:2 % employment this period
%             for j=1:2 % employment next period
%                 for n=1:2 % business cycle this period
%                     for m=1:2 % business cycle next period
%                         % capital choice next period from policy function
% %                         k2prime_bu=interpn(k,km,kprime(:,:,1,1),kprime,kmprime,'cubic'); 
% %                         k_next(i,:,j,n,m) = interp1(grid_k,k_guess(:,j,m),k_guess(:,i,n),'linear','extrap');
%                         k_next(i,:,:,j,n,m) = interpn(grid_k, km, k_guess(:,j,m), k_guess(:,i,n), kmprime, 'cubic');
% %                         k_next = interpn(grid_k, km, k_guess(:,:,j,m), k_guess, kmprime, 'cubic');
%                         % income matrix                         
%                         income = mat_income(K_guess,L(n),z(n));
%                      
%                         % consumption next period from budget constraint
%                         c_next(i,:,j,n,m) = max(1e-10,(1+r(K_guess,L(n),z(n))-delta)*k_guess(:,i,n) - [k_next(i,:,j,n,m)]' + income(:,m));
%                     end
%                 end
%             end
%         end

   % Bad aggregate state and unemployed idiosyncratic state 
   
     k_next_bu=interpn(grid_k,km,k_guess(:,:,1,1),k_guess,kmprime,'cubic'); 
        % finding the individual policy function k''=k(k',km') by interpolating
        % the previously found policy function k'=k(k,km) in new points (k',km')
%      cprime_bu=irate_b.*kprime+mu*(wage_b.*ones4)+(1-delta)*kprime-k_next_bu; 
                                                  % future consumption (c')
    
     c_next_bu=max(1e-10,(1+r_mat(km,L(1),z(1))-delta).*k_guess+mu*w_mat(km,L(1),z(1))-k_next_bu);                                            
     mu_next_bu=muc_inv(c_next_bu); % marginal utility of future consumption
   
   % Bad aggregate state and employed idiosyncratic state
     ones4 = ones(100,4,2,2);
     k_next_be=interpn(grid_k,km,k_guess(:,:,1,2),k_guess,kmprime,'cubic');
     c_next_be=max(1e-10,(1-delta+r_mat(km,L(1),z(1))).*k_guess+w_mat(km,L(1),z(1)).*((1*l_bar-mu*((U_b./(L_b))))*ones4) - k_next_be);
     mu_next_be=muc_inv(c_next_be);
 
   % Good aggregate state and unemployed idiosyncratic state
   
     k_next_gu=interpn(grid_k,km,k_guess(:,:,2,1),k_guess,kmprime,'cubic');
     c_next_gu=max(1e-10,(1+r_mat(km,L(2),z(2))-delta).*k_guess+mu*w_mat(km,L(2),z(2))-k_next_gu);
     mu_next_gu=muc_inv(c_next_gu);
   
   % Good aggregate state and employed idiosyncratic state
   
     k_next_ge=interpn(grid_k,km,k_guess(:,:,2,2),k_guess,kmprime,'cubic');
     c_next_ge=max(1e-10,(1-delta+r_mat(km,L(2),z(2))).*k_guess+w_mat(km,L(2),z(2)).*((1*l_bar-mu*((U_g./(L_g))))*ones4) - k_next_ge);
     mu_next_ge=muc_inv(c_next_ge);
        
        
%       calculate expected marginal utility of consumption next period
     Emuc_next =(mu_next_bu.*(1-delta+r_mat(km,L(1),z(1)))).*prob_bu+...
     (mu_next_be.*(1-delta+r_mat(km,L(1),z(1)))).*prob_be+...
     (mu_next_gu.*(1-delta+r_mat(km,L(1),z(1)))).*prob_gu+...
     (mu_next_ge.*(1-delta+r_mat(km,L(1),z(1)))).*prob_ge;
     
     c_current = muc_inv(beta*Emuc_next);
     k_current = NaN(100,4,2,2);
     raux = (1+r_mat(km,L(1),z(1))-delta);
     waux = w_mat(km,L(1),z(1));
     k_current(:,:,1,1) = raux(:,:,1,1).*grid_k+mu*waux(:,:,1,1)-c_current(:,:,1,1);
     k_current(:,:,1,2) = raux(:,:,1,2).*grid_k+waux(:,:,1,2).*((2*l_bar-mu*((U_b./(L_b))))*ones(100,4)) - c_current(:,:,1,2);
     k_current(:,:,2,1) = raux(:,:,2,1).*grid_k+mu*waux(:,:,2,1)-c_current(:,:,2,1);
     k_current(:,:,2,2) = raux(:,:,2,2).*grid_k+waux(:,:,2,2).*((2*l_bar-mu*((U_g./(L_g))))*ones(100,4)) - c_current(:,:,2,2);
     k_current=(k_current>=grid_k(1)).*(k_current<=grid_k(end)).*k_current+(k_current<grid_k(1))*k_min+(k_current>grid_k(end))*grid_k(end);
                       
       
        % apply borrowing constraint to get new policy function
%         k_new = max(k_min*K_guess,k_new); 
        d2=max(max(max(max(abs(k_new-k_guess)))));
%         d2_1 = norm(abs(k_new(:,1,1)-k_guess(:,1,1))./(1+abs(k_guess(:,1,1)))); % deviation between guess and new policy function
%         d2_2 = norm(abs(k_new(:,2,1)-k_guess(:,2,1))./(1+abs(k_guess(:,2,1))));
%         d2_3 = norm(abs(k_new(:,1,2)-k_guess(:,1,2))./(1+abs(k_guess(:,1,2))));
%         d2_4 = norm(abs(k_new(:,2,2)-k_guess(:,2,2))./(1+abs(k_guess(:,2,2))));
%         d2 = max([d2_1, d2_2, d2_3, d2_4]);
r_b = r_mat(km,L(1),z(1));
r_g = r_mat(km,L(2),z(2));
w_b = w_mat(km,L(1),z(1));
w_g = w_mat(km,L(2),z(2));
        save('some_results.mat', 'grid_k' ,'km','k_guess','kmprime', 'c_next_bu', 'mu_next_bu', 'k_next_bu', ...
            'c_next_be', 'mu_next_be', 'k_next_be', 'c_next_gu', 'mu_next_gu', 'k_next_gu', 'c_next_ge', 'mu_next_ge', 'k_next_ge', ...
            'Emuc_next', 'c_current', 'k_new', 'r_b', 'r_g', 'w_b', 'w_g');
        % update policy function
        k_guess = k_guess + 0.5*(k_new-k_guess);
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
    
%     K_demand = mean(mean(sim_k(ceil(T/2):end,:))); % average capital holdings over second half of sample
%     sim_L = mean(mean(id_shock(ceil(T/2):end,:)==2)); % average employment
%     
%     d1 = abs(K_demand-K_guess)./(1+K_guess); % deviation between guess for capital and demanded capital stock
%     
%     store.K_guess(iter)=K_guess;
%     store.K_next(iter)=K_demand;
% 
%     if K_demand>K_guess % update limits of bisection interval
%         K_lims(1) = K_guess;
%     else
%         K_lims(2) = K_guess;
%     end
%     
%     % update guess for aggregate capital stock
%     K_guess = (K_lims(1)+K_lims(2))/2;

% Time series for the ALM regression 

    ibad=0;           % count how many times the aggregate shock was bad
    igood=0;          % count how many times the aggregate shock was good
    xbad=0;  ybad=0;  % regression-variables for a bad state
    xgood=0; ygood=0; % regression-variables for a good state
    for i=ndiscard+1:T-1
        if agshock(i)==1
            ibad=ibad+1;
            xbad(ibad,1)=log(kmts(i));
            ybad(ibad,1)=log(kmts(i+1));
        else
            igood=igood+1;
            xgood(igood,1)=log(kmts(i));
            ygood(igood,1)=log(kmts(i+1));
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

%     if dif_B>(criter_B*100)
%         kcross=kcross1; % the new capital distribution  replaces the old one
%     end

    B=B1*update_B+B*(1-update_B); % update the vector of the ALM coefficients
    % according to the rule (9) in the paper
        
    disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K_guess),', K_demand: ',num2str(K_demand)])
end
toc