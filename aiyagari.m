clear
close all
%% set parameters
% consumers
beta = 0.99; % discount factor
sigma = 2; % risk aversion
mu = 0.3; % replacement rate of unemployed
k_min = 0.5; % borrowing constraint as share of capital stock

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
delta = 0.025; % depreciation rate
alpha = 0.33; % output elasticity of capital
z = 1; % productivity

% transition probabilities
L_target = 0.9;
PI_UE = 0.3; % chance of getting employed
PI_EU = PI_UE*(1-L_target)/L_target; % chance of getting unemployed


%% setup
% build transition matrix of agents
PI = [1-PI_UE,PI_UE;PI_EU,1-PI_EU]; % [UU,UE;EU;EE]
L = PI(1,2)/(PI(1,2)+PI(2,1)); % calculate steady state supply of labour
tau = mu*(1-L)/L; % tax rate

% define some useful functions
muc = @(c) c.^(-sigma); % marginal utility of consumption
muc_inv = @(muc) muc.^(-1/sigma); % inverse of marginal utility of consumption
w = @(K) z*(1-alpha)*K.^alpha*L^(-alpha); % wage
r = @(K) z*alpha*K.^(alpha-1)*L^(1-alpha); % rental rate of capital
K = @(r) L*(z*alpha/r).^(1/(1-alpha)); % capital stock necessary to yield return r
Y = @(K) z*K.^alpha*L^(1-alpha); % output
C = @(K) Y(K)-delta*K; % consumption

% define grid for individual capital on which to solve
K_rep = K(1/beta-1+delta);% capital stock of representative agent useful for comparison

grid_k_no = 100; % number of grid points for agents' capital holdings
grid_k = linspace(k_min*K_rep,3*K_rep,grid_k_no); % grid for agents' policy function

grid_dist_no = 1000;
grid_dist = linspace(grid_k(1),grid_k(end),grid_dist_no); % grid for distribution of agents - use finer grid

% useful matrices
mat_k = [grid_k',grid_k']; % replicate grid for unemployed and employed
mat_income = @(K) w(K)*repmat([mu,1-tau],grid_k_no,1); % matrix with income of each agent

%% solve for general equilibrium
HH_method = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
sim_method = 'histogram'; % 'simulation' for simulation, 'histogram' for histogram method
agg_method = 'bisection'; % 'bisection' to use bisection method, gradual updating otherwise

K_lims = [K_rep,grid_k(end)]; % initial limits for bisection method

% initial guesses
k_guess = mat_k; % policy function of agents (1st column unemployed, 2nd column employed)
if strcmp(agg_method,'bisection')
    K_guess = (K_lims(1)+K_lims(2))/2;
else
    K_guess = K_rep; % aggregate capital stock
end

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
        
        if strcmp(HH_method,'FP') % Fixed-point iteration     
            % calculate capital chosen next period, depending on current and
            % future employment status
            for i=1:2 % employment this period
                for j=1:2 % employment next period
                    % capital choice next period from policy function
                    k_next(i,:,j) = interp1(grid_k,k_guess(:,j),k_guess(:,i),'linear','extrap'); 
                end
                    % consumption next period from budget constraint
                    c_next(i,:,:) = max(1e-10,(1+r(K_guess)-delta)*[k_guess(:,i),k_guess(:,i)] - squeeze(k_next(i,:,:)) + mat_income(K_guess));
            end

            % calculate expected marginal utility of consumption next period
            Emuc_next(:,1) = PI(1,1)*muc(c_next(1,:,1))+PI(1,2)*muc(c_next(1,:,2));  % currently unemployed
            Emuc_next(:,2) = PI(2,1)*muc(c_next(2,:,1))+PI(2,2)*muc(c_next(2,:,2));  % currently employed       

            % calculate implied consumption this period from Euler equation
            c_current = muc_inv(beta*(1+r(K_guess)-delta)*Emuc_next);

            % calculate implied capital demand from budget constraint
            k_new = (1+r(K_guess)-delta)*mat_k + mat_income(K_guess) - c_current;
            
        elseif strcmp(HH_method,'FPend') % Fixed-point iteration with endogenous grid points method
            % calculate consumption next period, if agents are on the grid in
            % the next period
            c_next = max(1e-10,(1+r(K_guess)-delta)*mat_k-k_guess+mat_income(K_guess));

            % calculate expected marginal utility of consumption next period,
            % if future capital holdings are on the grid
            Emuc_next(:,1) = PI(1,1)*muc(c_next(:,1))+PI(1,2)*muc(c_next(:,2));  % currently unemployed
            Emuc_next(:,2) = PI(2,1)*muc(c_next(:,1))+PI(2,2)*muc(c_next(:,2));  % currently employed

            % calculate implied consumption this period from Euler equation
            c_current = muc_inv(beta*(1+r(K_guess)-delta)*Emuc_next);

            % calculate implied capital stock this period from budget constraint
            k_current = (c_current+mat_k-mat_income(K_guess))./(1+r(K_guess)-delta);

            % invert policy function to get k_next(k_current) on the grid for k
            for i=1:2 % currently unemployed and employed
                k_new(:,i) = interp1(k_current(:,i),grid_k,grid_k,'linear','extrap');
            end
        end

        % apply borrowing constraint to get new policy function
        k_new = max(k_min*K_guess,k_new); 
        
        d2 = norm(abs(k_new-k_guess)./(1+abs(k_guess))); % deviation between guess and new policy function
        
        % update policy function
        k_guess = k_guess + 0.5*(k_new-k_guess);
    end
    figure(1)
    plot(grid_k,log(k_guess./[grid_k',grid_k']))
    line([grid_k(1),grid_k(end)],[0,0])
    legend('unemployed','employed')
    xlabel('capital this period')
    ylabel('log(capital next period/capital this period)')
    title('policy functions')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Find distribution of agents
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if strcmp(sim_method ,'simulation')
        rng('default') % reset random number generator
        ind_no = 5000; % number of individuals simulated
        T = 5000; % number of periods simulated
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
        figure(4)
        subplot(1,2,1)
        plot(sum(sim_k,2)/ind_no) % plot aggregate capital demand
        xlabel('period')
        ylabel('capital demand')
        subplot(1,2,2)
        plot(sum(sim_e==2,2)/ind_no) % plot aggregate employment
        xlabel('period')
        ylabel('employment')    

        K_demand = mean(mean(sim_k(ceil(T/2):end,:))); % average capital holdings over second half of sample
        sim_L = mean(mean(sim_e(ceil(T/2):end,:)==2)); % average employment

    elseif strcmp(sim_method,'histogram')

        for i=1:2 % currently unemployed and employed
            k(:,i) = max(k_min*K_guess,min(grid_dist(end),interp1(grid_k,k_guess(:,i),grid_dist,'linear','extrap'))); % policy function on grid used for distribution

            temp = interp1(grid_dist,1:grid_dist_no,k(:,i),'linear','extrap'); % point in distribution to which weight is moved
            k_ind(:,i) = max(1,min(grid_dist_no-1,floor(temp))); % index of first grid point
            k_weight(:,i) = k_ind(:,i)+1-temp; % weight of first grid point

        end

        % build transition matrix of agents containing 4 blocks [UU,UE;EU;EE]
        % transition of capital for unemployed
        trans_U = sparse([1:grid_dist_no,1:grid_dist_no]',[k_ind(:,1);k_ind(:,1)+1],[k_weight(:,1);1-k_weight(:,1)],grid_dist_no,grid_dist_no);
        % transition of capital for employed
        trans_E = sparse([1:grid_dist_no,1:grid_dist_no]',[k_ind(:,2);k_ind(:,2)+1],[k_weight(:,2);1-k_weight(:,2)],grid_dist_no,grid_dist_no);

        transition = [PI(1,1)*trans_U,PI(1,2)*trans_U;PI(2,1)*trans_E,PI(2,2)*trans_E];
        warning('off','all')
        [eigvector,~] = eigs(transition',1,1); % find Eigenvector of transition matrix with Eigenvalue 1

        % alternative calculation of eigenvector by multiplying transition
        % matrix many times        
        %eigvector = (transition')^1000*ones(2*grid_dist_no,1);
        distribution = [eigvector(1:grid_dist_no),eigvector(grid_dist_no+(1:grid_dist_no))]/sum(eigvector);

        % alternative way to calculate distribution by solving linear
        % system of equations
        %distribution = reshape([(eye(2*grid_dist_no)-transition');ones(1,2*grid_dist_no,1)]\[zeros(2*grid_dist_no,1);1],grid_dist_no,2);

        K_demand = sum(sum(distribution.*k)); % aggregate capital demanded
        sim_L = sum(distribution(grid_dist_no+(1:grid_dist_no)));
    end
    
    d1 = abs(K_demand-K_guess)./(1+K_guess); % deviation between guess for capital and demanded capital stock
    
    store.K_guess(iter)=K_guess;
    store.K_next(iter)=K_demand;

    if K_demand>K_guess % update limits of bisection interval
        K_lims(1) = K_guess;
    else
        K_lims(2) = K_guess;
    end
    
    % update guess for aggregate capital stock
    if strcmp(agg_method,'bisection')
        K_guess = (K_lims(1)+K_lims(2))/2;
    else
        K_guess = K_guess + 0.05*(K_demand-K_guess);
    end
    disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K_guess),', K_demand: ',num2str(K_demand)])
end
toc

disp('_____________________________________________________________')
disp('Aggregate variables (log-deviation from representative agent)')
disp(['Capital:     ',num2str(K_guess),' (',num2str(log(K_guess/K_rep)),')'])
disp(['Output:      ',num2str(Y(K_guess)),' (',num2str(log(Y(K_guess)/Y(K_rep))),')'])
disp(['Consumption: ',num2str(C(K_guess)),' (',num2str(log(C(K_guess)/C(K_rep))),')'])


% Plot convergence of aggregate capital stock
figure(2)
plot([store.K_guess;store.K_next]')
legend('guess','demand')
xlabel('iteration')
ylabel('capital stock')

% Plot distribution of agents
figure(3)
if strcmp(sim_method,'simulation')
    temp = sim_k(ceil(T/2):end,:);
    hist(temp(:),100)
    legend('number of agents')
elseif strcmp(sim_method,'histogram')
    bar(grid_dist,distribution,'stacked')
    legend('unemployed','employed')
end
xlabel('capital holdings')
