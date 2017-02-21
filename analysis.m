%% Analysis
clear 
close all

%% Specify the type the solution method and simulation draws 
method.sim = 'histogram'; % specify solution method: 'histogram' or 'simulation'
par.ind_no = 100; % number of individuals simulated
par.T = 100; % number of periods simulated

%% Specify the parameter of interest and their values
method.analysis = 'Unemployment Benefit'; % specify parameter you want to analyze: 'Borrowing Constraint' for min savings, 'Unemployment Benefit', 'Chance to be Employed','Risk Aversion'
% par.vals = [0.15,0.25,0.35]; % indicate three parameter values for comparison
% par.vals2 = [0.4,0.25,0.2]; % for the second analysis, also change a second set of parameters

%% Plot grids
gridpar.start = 0.15;
gridpar.end = 0.35;
gridpar.no = 10;
mgrid.mu = linspace(gridpar.start,gridpar.end,gridpar.no);

gridpar.start2 = 0.4;
gridpar.end2 = 0.2;
mgrid.pi = linspace(gridpar.start2, gridpar.end2,gridpar.no);


%% Call the plotting function 
plotting_fct(par, method, gridpar, mgrid) 