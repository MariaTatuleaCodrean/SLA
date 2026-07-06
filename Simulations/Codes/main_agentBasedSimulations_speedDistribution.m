% Preamble
clear; addpath('Functions'); 

% Inputs
Nsim = 10;       % number of simulations
n0   = 1e4;      % number of bacteria
dx0  = 1e3;      % regular spacing between bacteria (1mm)

% Generate initial positions of bacteria (regular spacing)
X0 = (1:1:n0)*dx0;

% Stop simulations when we reach 45% active fraction
alpha_target = 0.3;
Nevent = round((1-alpha_target)*n0);

% Customize velocity distribution
for mf = [1/3 2/3 1]

    % Generate initial speed distribution
    muV = 9; sigV = mf*muV; 
    vel_pdf = generateDistribution(muV,sigV);

    % % Generate filename
    mf_str  = strrep(num2str(mf, '%.3f'), '.', 'p');
    outputfile = ['../Data/output_ABS_n_' num2str(n0) '_Nsim_' num2str(Nsim) '_sigV_per_muV_' mf_str '.mat'];
    save(outputfile,'dx0','mf','muV','sigV','alpha_target')

    % Simulations
    tic
    [x0,t0,tend,tsim,vel,len] = bactchainSimulation(n0,X0,vel_pdf,outputfile,Nsim,Nevent);
    toc
end