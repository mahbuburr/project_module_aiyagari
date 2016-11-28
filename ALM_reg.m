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

if dif_B>(criter_B*100)
    kcross=kcross1; % the new capital distribution  replaces the old one
end

B=B1*update_B+B*(1-update_B); % update the vector of the ALM coefficients 
                 % according to the rule (9) in the paper
iteration=iteration+1