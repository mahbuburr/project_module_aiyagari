function [ shock ] = prob_compare( probability, n )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

raux=rand;
if raux<=probability
    shock=n;
else
    shock=3-n;
end

end

