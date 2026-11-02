function [avg_vel, sem_vel, N] = ABSFitFunction_pooled(nvec,muV,sigV,n0,Nsim,tmin,tmax,hw,myseed)
%ABSFITFUNCTION Fit chain length vs. speed curve using predictions of agent
%based simulations (ABS).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   nvec        chain length (expressed as number of bacteria)
%   muV         mean of speed distribution
%   sigV        standard deviation of speed distribution
%   n0          number of bacteria
%   Nsim        number of simulations
%       * variables for spatio-temporal analysis *
%   tmin        start point 
%   tmax        end point 
%   hw          window half-width 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
%   avg_vel     mean chain speed for each chain length [um/s]
%   semv_vel    standard error of the mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Default parameters
if ~exist('seed','var') || isempty(myseed); myseed=2026; end

% Control parameters
dx0  = 1e3;         % regular spacing between bacteria (1mm)
X0 = (1:1:n0)*dx0;  % initial positions of bacteria (regular spacing)
xc = mean(X0);      % centre of initial array [um]

% Generate initial speed distribution
vel_pdf = generateDistribution(muV,sigV);

% Simulations
obs_vel = cell(1,Nsim);
obs_len = cell(1,Nsim);
tic
for jj = 1:Nsim
    % Independent, reproducible random stream for each lane (Claude Opus 5.5 suggestion)
    % To get a fresh, independent set of simulations (for example, to check 
    % that a result isn’t a fluke of one seed), change 2026 to another number.
    s = RandStream('Threefry','Seed',myseed);
    s.Substream = jj;
    RandStream.setGlobalStream(s);
    
    % Run simulation
    [x0,t0,tend,~,vel,len,way] = bactchainSimulation_pooled(n0,X0,vel_pdf,tmax,1);

    % Determine which agents were observed in the spatio-temporal window
    valid = ~isnan(t0);                  % slots that were actually used

    tend(isnan(tend)) = Inf;             % chains still alive at tmax
    ta = max(t0,tmin);                   % portion of life inside [tmin,tmax]
    tb = min(tend,tmax);
    inTime  = ta <= tb;                  % inside time window?

    xa = x0 + (ta-t0).*vel.*way;         % position at ta
    xb = x0 + (tb-t0).*vel.*way;         % position at tb
    inSpace = (min(xa,xb) <= xc+hw) ...  
        & (max(xa,xb) >= xc-hw);         % inside space window?

    obs = valid & inTime & inSpace;      % observe valid chains in space & time

    obs_vel{jj} = vel(obs);
    obs_len{jj} = len(obs);

end
toc

% Limit analysis to chains of length...
vel_all = [obs_vel{:}];
len_all = [obs_len{:}];

% Define variables
avg_vel = nan(size(nvec));
sem_vel = nan(size(nvec));

% Analyze results
for kk = 1:numel(nvec)
    v = vel_all(len_all==nvec(kk));
    avg_vel(kk) = mean(v);
    sem_vel(kk) = std(v)/sqrt(numel(v));
end

% Report number of detected agents
N = numel(vel_all);

end