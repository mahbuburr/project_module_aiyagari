clear all
close all
clc

parameters; % load parameters

[id_shock,ag_shock]  = generate_shocks(prob,T,ind_no,U_b); % generate shocks