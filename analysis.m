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
% Create grid of unemployment benefits for the analysis
gridpar.start = 0.01;
gridpar.end = 0.6;
gridpar.no = 20;
mgrid.mu = linspace(gridpar.start, gridpar.end, gridpar.no);
mgrid.mu(7) = mgrid.mu(6); % the original point does not converge
mgrid.mu(6) = 0.15; % the original point does not converge
mgrid.mu(8) = 0.18; % the original point does not converge

% matching grid of the transition probabilities
gridpar.start2 = 0.418006431;
gridpar.end2 = 0.351351351;
mgrid.pi = [0.418006431, 0.413867753, 0.409823146, 0.405844156, 0.401954115,0.4, 0.398125746, 0.396341463, 0.38708909, 0.383537395, 0.380061395, 0.376636922, 0.373284328, 0.369980363, 0.366744717, 0.363555009, 0.360430298, 0.35734902, 0.354329636, 0.351351351];


%% Call the plotting function 
plotting_fct(par, method, gridpar, mgrid) 