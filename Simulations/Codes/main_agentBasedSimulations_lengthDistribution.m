% Preamble
clear; addpath('Functions') % load user defined functions

% Inputs
Nsim = 10;       % number of simulations
n0   = 1e4;      % number of bacteria 
dx0  = 1e3;      % regular spacing between bacteria (1mm)

% Generate initial positions of bacteria (regular spacing)
X0 = (1:1:n0)*dx0;

% Generate initial speed distribution
muV  = 9.0; sigV = 6.0;
vel_pdf = generateDistribution(muV,sigV);

% Stop simulations when we reach 45% active fraction
alpha_target = 0.25;
Nevent = round((1-alpha_target)*n0);

% Customize length distribution
for sigL = [1/3 2/3 1]

    % Generate bacterial length distribution
    muL = 1;
    len_pdf = generateDistribution(muL,sigL);
    % len_pdf = makedist('Uniform','lower',2/3,'upper',4/3);

    % Generate filename
    sigL_str  = strrep(num2str(sigL, '%.3f'), '.', 'p');
    outputfile = ['../Data/output_ABS_n_' num2str(n0) '_Nsim_' num2str(Nsim) '_sigL_per_muL_' sigL_str 'Lognormal.mat'];
    save(outputfile,'dx0','sigL','muV','sigV','alpha_target')

    % Simulations
    tic
    [x0,t0,tend,tsim,vel,len,way,active,alpha,nbact] = bactchainSimulation_lengthdistrib(n0,X0,vel_pdf,len_pdf,outputfile,Nsim,Nevent);
    toc

    % Analyze results
    tplot = tsim(1,end);    % select simulation end for plotting
    nmax = 8;               % limit analysis to chains of length...

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
            I = find((round(len(jj,:))==m) & active(jj,:));

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
end