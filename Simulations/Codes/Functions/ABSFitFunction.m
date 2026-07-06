function [norm_v] = ABSFitFunction(nvec,muV,sigV)
%ABSFITFUNCTION Fit chain length vs. speed curve using predictions of agent
%based simulations (ABS).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   nvec        chain length (expressed as number of bacteria)
%   muV         mean of speed distribution
%   sigV        variance of speed distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
%   norm_v      normalized velocity for each chain length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Initialize output variables (same size as nvec)
norm_v  = zeros(size(nvec));

% Control parameters
Nsim = 10;          % number of simulations
n0   = 1e4;         % number of bacteria 
dx0  = 1e3;         % regular spacing between bacteria (1mm)
X0 = (1:1:n0)*dx0;  % initial positions of bacteria (regular spacing)

% Generate initial speed distribution
vel_pdf = generateDistribution(muV,sigV);

% Stop simulations when we reach 45% active fraction
alpha_target = 0.45;
Nevent = round((1-alpha_target)*n0);

% Generate filename & file
outputfile = 'autosave.mat';
save(outputfile,'Nsim','n0','dx0','muV','sigV')

% Simulations
[~,t0,tend,tsim,vel,len] = bactchainSimulation(n0,X0,vel_pdf,outputfile,Nsim,Nevent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include chains active at final step
tend(isnan(tend)) = Inf;   

% For each simulation
for jj=1:Nsim
    % Select last time point of simulation (when we reached alpha_target)
    tstar = tsim(1,end);

    % Calculate average velocity of single bacteria
    V1 = mean(vel(jj,...
        (t0(jj,:)<tstar) & (tend(jj,:)>tstar) & (len(jj,:)==1)));

    for ii = 1:length(nvec)
        % Calculate average velocity of chains of length nvec(ii)
        avg_vel = mean(vel(jj,...
            (t0(jj,:)<tstar) & (tend(jj,:)>tstar) & (len(jj,:)==nvec(ii))));

        % Calculate normalized speed
        norm_v(ii) = norm_v(ii) + avg_vel; %/V1;
    end  
end

% Average over all simulations
norm_v = norm_v/Nsim;

% Display current CPU time
currentCPUTime = cputime;
fprintf('Current CPU time: %.6f seconds\n', currentCPUTime);
end