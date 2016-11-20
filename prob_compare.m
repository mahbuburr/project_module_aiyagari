function [ shock ] = prob_compare( p, n )
% prob_compare(p, n) draws a random number which compares to probability, p,
% and returns n if the number is less or equal to p. Otherwise it returns
% 3-n. 
% 
% n could take values 1 or 2.
% p is between -1 and 1.

if rand<=p
    shock=n;
else
    shock=3-n;
end

end

