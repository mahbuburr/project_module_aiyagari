clear all
close all
tic

parfor i = 1:100000000
   B(i) = 10^i;
end
toc