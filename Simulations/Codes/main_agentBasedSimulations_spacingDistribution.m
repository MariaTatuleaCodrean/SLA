% Preamble
clear; addpath('Functions'); 

% Inputs
Nsim = 10;       % number of simulations
n0   = 1e4;      % number of bacteria         
dx0  = 1e3;      % regular spacing between bacteria (1mm)

% Generate initial speed distribution
muV = 9; sigV = 6; vel_pdf = generateDistribution(muV,sigV);

% Stop simulations when we reach 45% active fraction
alpha_target = 0.3;
Nevent = round((1-alpha_target)*n0);

% Customize spacing distribution
for mf = [0 0.3 0.6 0.9] % size of displacement interval

    % Generate initial positions of bacteria (irregular spacing)
    dX0 = unifrnd(-.5*mf*dx0,.5*mf*dx0,[1 n0]);
    X0  = (1:1:n0)*dx0 + dX0;

    % Generate filename
    mf_str  = strrep(num2str(mf, '%.3f'), '.', 'p');
    outputfile = ['../Data/output_ABS_n_' num2str(n0) '_Nsim_' num2str(Nsim) '_uniform_' mf_str 'd0_interval.mat'];
    save(outputfile,'dx0','mf','dX0','muV','sigV','alpha_target')

    % Simulations
    tic
    [x0,t0,tend,tsim,vel,len] = bactchainSimulation(n0,X0,vel_pdf,outputfile,Nsim,Nevent);
    toc
end