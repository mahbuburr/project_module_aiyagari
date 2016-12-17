function [ prob, prob_aux ] = probability_matrix( U_b, PI_UE_b, U_g, PI_UE_g, B_p, PI_bg, grid_k_no, grid_K_no, ag_states_no, id_states_no)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

L_b=(1-U_b);    % employment rate in a bad aggregate state
PI_EU_b = PI_UE_b*U_b/L_b; % chance of getting unemployed in bad state
PI_b = [1-PI_UE_b,PI_UE_b;PI_EU_b,1-PI_EU_b]; % [UU,UE;EU;EE]

L_g=(1-U_g);    % employment rate in a good aggregate state
PI_EU_g = PI_UE_g*U_g/L_g; % chance of getting unemployed in good state
PI_g = [1-PI_UE_g,PI_UE_g;PI_EU_g,1-PI_EU_g]; % [UU,UE;EU;EE]

G_p = 1 - B_p; % good states proportion
PI_gb = PI_bg*B_p/G_p; % chance of observing good state when at bad state
PI = [1-PI_bg, PI_bg;PI_gb,1-PI_gb];

prob = NaN(4,4);
prob(1:2,1:2) = PI_b.*PI(1,1);
prob(3:4,3:4) = PI_g.*PI(2,2);
prob(3,1) = PI(2,1)*1.25*prob(1,1)/PI(1,1);
prob(3,2) = PI(2,1) - prob(3,1);
prob(4,1) = (U_b*PI(2,1)-U_g*prob(3,1))/L_g;
prob(4,2) = PI(2,1) - prob(4,1);
prob(1,3) = PI(1,2)*0.75*prob(3,3)/PI(2,2);
prob(1,4) = PI(1,2) - prob(1,3);
prob(2,3) = (U_g*PI(1,2)-U_b*prob(1,3))/L_b;
prob(2,4) = PI(1,2) - prob(2,3);


% Augumented probability matrix 
% n = 1 bad agg. state and unemployed in next period
% n = 2 bad agg. state and employed in next period
% n = 3 good agg. state and unemployed in next period
% n = 4 good agg. state and employed in next period
prob_aux = zeros(grid_k_no,grid_K_no,ag_states_no,id_states_no, 4);
for n = 1:4
    m=1;
    for i = 1:2
        for j = 1:2
            prob_aux(:,:,i,j,n) = repmat(prob(m,n),grid_k_no,grid_K_no);
            m=m+1;
        end
    end
end

end
