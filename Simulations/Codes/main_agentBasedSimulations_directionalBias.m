% Preamble
clear; addpath('Functions'); 

% Inputs
Nsim = 1;       % number of simulations
n0   = 1e4;      % number of bacteria
dx0  = 1e3;      % regular spacing between bacteria (1mm)

% Generate initial positions of bacteria (regular spacing)
X0 = (1:1:n0)*dx0;

% Generate initial speed distribution
muV  = 9.0; sigV = 3.0;
vel_pdf = generateDistribution(muV,sigV);

% Stop simulations when we reach 45% active fraction
alpha_target = 0.05;
Nevent = round((1-alpha_target)*n0);

% Customize directional bias
for Pright = [.5 .75 1]

    % Generate filename
    Pright_str  = strrep(num2str(Pright, '%.2f'), '.', 'p');
    outputfile = ['../Data/output_ABS_n_' num2str(n0) '_Nsim_' num2str(Nsim) '_Pright_' Pright_str 'muV9_sigV3.mat'];
    save(outputfile,'dx0','Pright','muV','sigV','alpha_target')

    % Simulations
    tic
    [x0,t0,tend,tsim,vel,len] = bactchainSimulation(n0,X0,vel_pdf,outputfile,Nsim,Nevent,Pright);
    toc
end