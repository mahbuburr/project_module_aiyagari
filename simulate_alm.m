function [ output_args ] = simulate_alm( K_demand, ag_shock, B, T )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

K_alm=zeros(T,1);  % represents aggregate capital computed from the ALM
K_alm(1)=K_demand(1);  % in the first period km computed from the ALM (kmalm) 
                   % is equal km computed from the cross-sectional capital 
                   % distribution (kmts)
                                   
for t=1:T-1       % compute kmalm for t=2:T
   if ag_shock(t)==1
      K_alm(t+1)=exp(B(1)+B(2)*log(K_alm(t)));
   else
      K_alm(t+1)=exp(B(3)+B(4)*log(K_alm(t)));
   end
   
end

output_args = mean(K_alm(ceil(T/2):end));
end

