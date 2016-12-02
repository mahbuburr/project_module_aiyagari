function [ prob ] = probability_matrix( U_b, PI_UE_b, U_g, PU_UE_g, B_p, PI_bg)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

U_b=0.1;        % unemployment rate in a bad aggregate state
L_b=(1-U_b);    % employment rate in a bad aggregate state
PI_UE_b = 0.4;   % chance of getting employed in bad state
PI_EU_b = PI_UE_b*U_b/L_b; % chance of getting unemployed in bad state
PI_b = [1-PI_UE_b,PI_UE_b;PI_EU_b,1-PI_EU_b]; % [UU,UE;EU;EE]

U_g=0.04;       % unemployment rate in a good aggregate state
L_g=(1-U_g);    % employment rate in a good aggregate state
PI_UE_g = 2/3;  % chance of getting employed in good state
PI_EU_g = PI_UE_g*U_g/L_g; % chance of getting unemployed in good state
PI_g = [1-PI_UE_g,PI_UE_g;PI_EU_g,1-PI_EU_g]; % [UU,UE;EU;EE]

L = [L_b; L_g]; % vector of states for labour
l_bar=1/L_b;    % used for simplification

B_p = 0.5; % bad states proportion
G_p = 1 - B_p; % good states proportion
PI_bg = 0.125; % chance of observing a bad state when at good state
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

end
