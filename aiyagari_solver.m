function [ k, c, K, sim, store] = aiyagari_solver( par, func, method)
% AIYAGARI MODEL: Heterogeneous agents model due to idiosyncratic labour
% shocks. Agents self-sinsure against unemploment by building capital
% stock.
% Input variables:
%   par = (calibrated) Parameters that decribe the economy.
%   grid = Grid to calculate policy functions and finer grid for simulation
%       economy.
%   K = Aggregate capital. Include starting guess and solution for
%       representative consumer.
%   k = Individual capital. Also include a starting guess.
%   func = Helpful functions to automate certain calculations.
%   method = Describe the method of iteration/ simutation.
%   mat = Grid and income for unemployed and employed in one matrice each.
% 
% Output variables:

% define grid for individual capital on which to solve
K.rep = func.K(1/par.beta-1+par.delta);% capital stock of representative agent useful for comparison

grid.k_no = 100; % number of grid points for agents' capital holdings
grid.k = linspace(par.k_min*K.rep,3*K.rep,grid.k_no); % grid for agents' policy function

grid.dist_no = 1000;
grid.dist = linspace(grid.k(1),grid.k(end),grid.dist_no); % grid for distribution of agents - use finer grid

% useful matrices
mat.k = [grid.k',grid.k']; % replicate grid for unemployed and employed
mat.income = @(K) func.w(K)*repmat([par.mu,1-par.tau],grid.k_no,1); % matrix with income of each agent

K.lims = [K.rep,grid.k(end)]; % initial limits for bisection method

% initial guesses
k.guess = mat.k; % policy function of agents (1st column unemployed, 2nd column employed)
if strcmp(method.agg,'bisection')
    K.guess = (K.lims(1)+K.lims(2))/2;
else
    K.guess = K.rep; % aggregate capital stock
end

d1 = 1;
iter = 0;
while d1>1e-6 && iter<50 % loop for aggregate problem
    iter = iter+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Solve household problem given prices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d2 = 1;
    while d2>1e-10 % loop for household problem
        
        if strcmp(method.HH,'FP') % Fixed-point iteration     
            % calculate capital chosen next period, depending on current and
            % future employment status
            for i=1:2 % employment this period
                for j=1:2 % employment next period
                    % capital choice next period from policy function
                    k.next(i,:,j) = interp1(grid.k,k.guess(:,j),k.guess(:,i),'linear','extrap'); 
                end
                    % consumption next period from budget constraint
                    c.next(i,:,:) = max(1e-10,(1+func.r(K.guess)-par.delta)*[k.guess(:,i),k.guess(:,i)] - squeeze(k.next(i,:,:)) + mat_income(K_guess));
            end

            % calculate expected marginal utility of consumption next period
            Emuc_next(:,1) = par.PI(1,1)*func.muc(c.next(1,:,1))+par.PI(1,2)*func.muc(c.next(1,:,2));  % currently unemployed
            Emuc_next(:,2) = par.PI(2,1)*func.muc(c.next(2,:,1))+par.PI(2,2)*func.muc(c.next(2,:,2));  % currently employed       

            % calculate implied consumption this period from Euler equation
            c.current = func.muc_inv(par.beta*(1+func.r(K.guess)-par.delta)*Emuc_next);

            % calculate implied capital demand from budget constraint
            k.new = (1+func.r(K.guess)-par.delta)*mat.k + mat.income(K.guess) - c.current;
            
        elseif strcmp(method.HH,'FPend') % Fixed-point iteration with endogenous grid points method
            % calculate consumption next period, if agents are on the grid in
            % the next period
            c.next = max(1e-10,(1+func.r(K.guess)-par.delta)*mat.k-k.guess+mat.income(K.guess));

            % calculate expected marginal utility of consumption next period,
            % if future capital holdings are on the grid
            Emuc_next(:,1) = par.PI(1,1)*func.muc(c.next(:,1))+par.PI(1,2)*func.muc(c.next(:,2));  % currently unemployed
            Emuc_next(:,2) = par.PI(2,1)*func.muc(c.next(:,1))+par.PI(2,2)*func.muc(c.next(:,2));  % currently employed

            % calculate implied consumption this period from Euler equation
            c.current = func.muc_inv(par.beta*(1+func.r(K.guess)-par.delta)*Emuc_next);

            % calculate implied capital stock this period from budget constraint
            k.current = (c.current+mat.k-mat.income(K.guess))./(1+func.r(K.guess)-par.delta);

            % invert policy function to get k_next(k_current) on the grid for k
            for i=1:2 % currently unemployed and employed
                k.new(:,i) = interp1(k.current(:,i),grid.k,grid.k,'linear','extrap');
            end
        end

        % apply borrowing constraint to get new policy function
        k.new = max(par.k_min*K.guess,k.new); 
        
        d2 = norm(abs(k.new-k.guess)./(1+abs(k.guess))); % deviation between guess and new policy function
        
        % update policy function
        k.guess = k.guess + 0.5*(k.new-k.guess);
    end
    figure(1)
    plot(grid.k,log(k.guess./[grid.k',grid.k']))
    line([grid.k(1),grid.k(end)],[0,0])
    legend('unemployed','employed')
    xlabel('capital this period')
    ylabel('log(capital next period/capital this period)')
    title('policy functions')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    Find distribution of agents
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if strcmp(method.sim ,'simulation')
        rng('default') % reset random number generator
        ind_no = 5000; % number of individuals simulated
        T = 5000; % number of periods simulated
        sim.k = zeros(T,ind_no); % simulated values of capital stock
        sim.e = ones(T,ind_no); % simulated employment status
        sim.shock = rand(T,ind_no); % shocks for employment transition

        sim.e(1,1:round(par.L*ind_no))=2; % initial individuals that are employed
        sim.k(1,:) = K.guess; % initial capital holdings

        for t=2:T
            sim.e(t,sim.e(t-1,:)==1) = 1+(sim.shock(t,sim.e(t-1,:)==1)<=par.PI(1,2)); % new employment status of previously unemployed
            sim.e(t,sim.e(t-1,:)==2) = 1+(sim.shock(t,sim.e(t-1,:)==2)<=par.PI(2,2)); % new employment status of previously employed    
            sim.k(t,sim.e(t,:)==1) = interp1(grid.k,k.guess(:,1),sim.k(t-1,sim.e(t,:)==1),'linear','extrap'); % capital demand of currently unemployed
            sim.k(t,sim.e(t,:)==2) = interp1(grid.k,k.guess(:,2),sim.k(t-1,sim.e(t,:)==2),'linear','extrap'); % capital demand of currently employed
        end
        figure(4)
        subplot(1,2,1)
        plot(sum(sim.k,2)/ind_no) % plot aggregate capital demand
        xlabel('period')
        ylabel('capital demand')
        subplot(1,2,2)
        plot(sum(sim.e==2,2)/ind_no) % plot aggregate employment
        xlabel('period')
        ylabel('employment')    

        K.demand = mean(mean(sim.k(ceil(T/2):end,:))); % average capital holdings over second half of sample
        sim.L = mean(mean(sim.e(ceil(T/2):end,:)==2)); % average employment

    elseif strcmp(method.sim,'histogram')

        for i=1:2 % currently unemployed and employed
            k.k(:,i) = max(par.k_min*K.guess,min(grid.dist(end),interp1(grid.k,k.guess(:,i),grid.dist,'linear','extrap'))); % policy function on grid used for distribution

            temp = interp1(grid.dist,1:grid.dist_no,k.k(:,i),'linear','extrap'); % point in distribution to which weight is moved
            k.ind(:,i) = max(1,min(grid.dist_no-1,floor(temp))); % index of first grid point
            k.weight(:,i) = k.ind(:,i)+1-temp; % weight of first grid point

        end

        % build transition matrix of agents containing 4 blocks [UU,UE;EU;EE]
        % transition of capital for unemployed
        trans_U = sparse([1:grid.dist_no,1:grid.dist_no]',[k.ind(:,1);k.ind(:,1)+1],[k.weight(:,1);1-k.weight(:,1)],grid.dist_no,grid.dist_no);
        % transition of capital for employed
        trans_E = sparse([1:grid.dist_no,1:grid.dist_no]',[k.ind(:,2);k.ind(:,2)+1],[k.weight(:,2);1-k.weight(:,2)],grid.dist_no,grid.dist_no);

        transition = [par.PI(1,1)*trans_U,par.PI(1,2)*trans_U;par.PI(2,1)*trans_E,par.PI(2,2)*trans_E];
        warning('off','all')
        [eigvector,~] = eigs(transition',1,1); % find Eigenvector of transition matrix with Eigenvalue 1

        % alternative calculation of eigenvector by multiplying transition
        % matrix many times        
        %eigvector = (transition')^1000*ones(2*grid_dist_no,1);
        store.distribution = [eigvector(1:grid.dist_no),eigvector(grid.dist_no+(1:grid.dist_no))]/sum(eigvector);

        % alternative way to calculate distribution by solving linear
        % system of equations
        %distribution = reshape([(eye(2*grid_dist_no)-transition');ones(1,2*grid_dist_no,1)]\[zeros(2*grid_dist_no,1);1],grid_dist_no,2);
        K.demand = sum(sum(store.distribution.*k.k)); % aggregate capital demanded
        sim.L = sum(store.distribution(grid.dist_no+(1:grid.dist_no)));
    end
    
    d1 = abs(K.demand-K.guess)./(1+K.guess); % deviation between guess for capital and demanded capital stock
    
    store.K_guess(iter)=K.guess;
    store.K_next(iter)=K.demand;

    if K.demand>K.guess % update limits of bisection interval
        K.lims(1) = K.guess;
    else
        K.lims(2) = K.guess;
    end
    
    % update guess for aggregate capital stock
    if strcmp(method.agg,'bisection')
        K.guess = (K.lims(1)+K.lims(2))/2;
    else
        K.guess = K.guess + 0.05*(K.demand-K.guess);
    end
    disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K.guess),', K_demand: ',num2str(K.demand)])
end


disp('_____________________________________________________________')
disp('Aggregate variables (log-deviation from representative agent)')
disp(['Capital:     ',num2str(K.guess),' (',num2str(log(K.guess/K.rep)),')'])
disp(['Output:      ',num2str(func.Y(K.guess)),' (',num2str(log(func.Y(K.guess)/func.Y(K.rep))),')'])
disp(['Consumption: ',num2str(func.C(K.guess)),' (',num2str(log(func.C(K.guess)/func.C(K.rep))),')'])


% Plot convergence of aggregate capital stock
figure(2)
plot([store.K_guess;store.K_next]')
legend('guess','demand')
xlabel('iteration')
ylabel('capital stock')

% Plot distribution of agents
figure(3)
if strcmp(method.sim,'simulation')
    temp = sim.k(ceil(T/2):end,:);
    hist(temp(:),100)
    legend('number of agents')
elseif strcmp(method.sim,'histogram')
    bar(grid.dist,store.distribution,'stacked')
    legend('unemployed','employed')
end
xlabel('capital holdings')

end

