function [ k_guess ] = solve_individual_problem( mu, k_guess, K_guess )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Solve household problem given prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixed_parameters;
tau = mu*(1-L)/L; % tax rate
mat_income = @(K) w(K)*repmat([mu,1-tau],grid_k_no,1); % matrix with income of each agent
d2 = 1;
while d2>1e-10 % loop for household problem
    
    
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
    c_current = muc_inv(betta*(1+r(K_guess)-delta)*Emuc_next);
    
    % calculate implied capital demand from budget constraint
    k_new = (1+r(K_guess)-delta)*mat_k + mat_income(K_guess) - c_current;
        
    % apply borrowing constraint to get new policy function
    k_new = max(k_min*K_guess,k_new);
    
    d2 = norm(abs(k_new-k_guess)./(1+abs(k_guess))); % deviation between guess and new policy function
    
    % update policy function
    k_guess = k_guess + 0.5*(k_new-k_guess);
end
    