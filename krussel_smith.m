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
    while d2>1e-10 % loop for household problem
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


        
        
%       calculate expected marginal utility of consumption next period
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
                   % calculate implied consumption this period from Euler equation
                   c_current(:,i,n) = muc_inv(beta*(1+r(K_guess,L(m),z(m))-delta)*Emuc_next(:,i,n));
                   % income matrix
                   income = mat_income(K_guess,L(n),z(n));
                   
                   % calculate implied capital demand from budget constraint
                   k_new(:,i,n) = (1+r(K_guess,L(n),z(n))-delta)*grid_k' + income(:,n) - c_current(:,i,n);
               end
           end
       end
                       
       
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