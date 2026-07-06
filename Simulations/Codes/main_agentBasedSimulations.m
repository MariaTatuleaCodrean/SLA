% Preamble
clear; addpath('Functions'); 

% Inputs
Nsim = 10;          % number of simulations
n0   = 1e4;         % number of bacteria 
dx0  = 1e3;         % regular spacing between bacteria (1mm)
X0 = (1:1:n0)*dx0;  % initial positions of bacteria (regular spacing)

% User input
choice = 'fig1e_parameters';

% Generate initial speed distribution
switch true
    case strcmp(choice,'fig1e_parameters')
        load('../Data/modelfit_30March2026.mat','muV','sigV')
        vel_pdf = generateDistribution(muV,sigV);

    case strcmp(choice,'fig3bc_parameters')
        muV = 9; sigV = 6;
        vel_pdf = generateDistribution(muV,sigV);
end

% Stop simulations when we reach a given active fraction
switch true
    case strcmp(choice,'fig1e_parameters')
        alpha_target = 0.45;
    case strcmp(choice,'fig3bc_parameters')
        alpha_target = 0.25; 
end

% Compute number of merging events
Nevent = round((1-alpha_target)*n0);

% Generate filename
muV_str  = strrep(num2str(muV, '%.3f'), '.', 'pt');
sigV_str  = strrep(num2str(sigV, '%.3f'), '.', 'pt');
outputfile = ['../Data/output_ABS_n_' num2str(n0) '_Nsim_' num2str(Nsim) '_muV_' muV_str '_sigV_' sigV_str '.mat']; 
save(outputfile,'Nsim','n0','dx0','X0','muV','sigV','vel_pdf','alpha_target','Nevent')

% Simulations
tic
[x0,t0,tend,tsim,vel,len,way,active,alpha] = bactchainSimulation(n0,X0,vel_pdf,outputfile,Nsim,Nevent);
toc

%% Analyze results (OPTIONAL)
tplot = mean(tsim(:,end));  % select simulation end for plotting     
nmax  = 8;                  % limit analysis to chains of length...

% Define variables
avg_vel = nan(nmax,Nsim);
std_vel = nan(nmax,Nsim);
sem_vel = nan(nmax,Nsim);
sim_vel = cell(1,Nsim);
sim_len = cell(1,Nsim);

% Determine evolution of chain length and speed
for jj=1:Nsim 
    % Identify all active chains at end of simulation
    I = find(active(jj,:));
    sim_vel{jj} = vel(jj,I); 
    sim_len{jj} = len(jj,I);

    % For each chain length
    for m = 1:nmax 
        % Identify active chains of length m 
        I = find((len(jj,:)==m) & active(jj,:));

        % Calculate average velocity
        avg_vel(m,jj) = mean(vel(jj,I));

        % Calculate standard deviation
        w = 0; %  the standard deviation is normalized by N-1
        std_vel(m,jj) = std(vel(jj,I),w);

        % Calculate standard error of the mean
        % Recall that s.e.m. = std(x)/sqrt(length(x));
        sem_vel(m,jj) = std_vel(m,jj)/sqrt(length(I));
    end
end

save(outputfile,'tplot','nmax','sim_vel','sim_len','avg_vel','sem_vel','std_vel','-append')