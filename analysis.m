%% Analysis

clear 
close all

%% Specify the type the solution method and simulation draws 
method.sim = 'histogram'; % specify solution method: 'histogram' or 'simulation'
par.ind_no = 100; % number of individuals simulated
par.T = 100; % number of periods simulated

%% Specify the parameter of interest and their values
method.analysis = 'Chance to be Employed'; % specify parameter you want to analyze: 'Borrowing Constraint' for min savings, 'Unemployment Benefit', 'Chance to be Employed','Risk Aversion'
par.vals = [0.1,0.3,0.9]; %indicate three parameter values for comparison

%% Call the plotting function 
plotting_fct(par, method) 