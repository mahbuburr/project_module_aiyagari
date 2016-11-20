function [id_shock,ag_shock]  = generate_shocks(prob,T,N,U_b)

disp('Generating shocks');

id_shock=zeros(T,N); % matrix of idiosyncratic shocks 
ag_shock=zeros(T,1); % vector of aggregate shocks
prob_block = mat2cell(prob, [2 2], [2 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition probababilities between the aggregate states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prob_ag(i,j) is the probability of tomorrow's agg. shock (i=1,2) given 
% today's agg. shock (j=1,2)

prob_ag=NaN(2,2);
prob_ag(1,1)=prob(1,1)+prob(1,2); prob_ag(2,1)=1-prob_ag(1,1);  
prob_ag(2,2)=prob(3,3)+prob(3,4); prob_ag(1,2)=1-prob_ag(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition probabilities of an idiosyncratic shock epsilon' given that 
% aggregate shock s' is realized;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P(x,y,z,w) is the probability of idiosyncratic shock epsilon'=w conditional 
% on aggregate shocks s'=y, s=x and idiosyncratic shock epsilon=z

P = NaN(2,2,2,2);
for i = 1:2
    for j = 1:2
        temp = diag(prob_block{i,j})./prob_ag(j,i);
        P(j,i,1,1) = temp(1);
        P(j,i,1,2) = 1-temp(1);
        P(j,i,2,2) = temp(2);
        P(j,i,2,1) = 1 - temp(2);
    end
end

%__________________________________________________________________________
%
% Generation of the aggregate shocks 
%__________________________________________________________________________

ag_shock(1)=1; % assume that initially (t=1), the economy is in a bad state 

% To generate shocks in subsequent periods, we draw random numbers. If a 
% random number drawn is <= probability of a bad shock conditional on 
% agshock(t-1)), then set agshock(t)=1; otherwise set agshock(t)=2

for t=2:T
   ag_shock(t)=prob_compare(prob_ag(1,ag_shock(t-1)),1);
end

%__________________________________________________________________________
%
% Generation of the idiosyncratic shocks for all agents in the first period
%__________________________________________________________________________

for i=1:N
      id_shock(1,i)=prob_compare(U_b, 1);
end
%__________________________________________________________________________
%
% Generation of the idiosyncratic shocks for all agents starting from the 
% second period
%__________________________________________________________________________

% P(s,s',e,e') is the probability of idiosyncratic shock epsilon'=w conditional 
% on aggregate shocks s'=y, s=x and idiosyncratic shock epsilon=z

for t = 2:T
    for q = 1:2
        for p = 1:2
            if ag_shock(t-1)==q && ag_shock(t)==p
                for i=1:N
                    for j = 1:2
                        if id_shock(t-1,i)==j
                            id_shock(t,i) = prob_compare(P(q,p,j,j),j);
                        end
                    end
                end
            end
        end
    end
end
                         
end