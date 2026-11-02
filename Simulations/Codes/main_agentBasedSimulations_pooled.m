% Preamble
clear; addpath('Functions'); 

% Inputs
Nsim = 1000;          % number of independent single-lane simulations
n0   = 20;          % number of bacteria 
dx0  = 1e3;         % regular spacing between bacteria (1mm)
X0 = (1:1:n0)*dx0;  % initial positions of bacteria (regular spacing)

% Generate initial speed distribution - FIG. 1E parameters
load('../Data/modelfit_29October2026.mat','muV','sigV')
vel_pdf = generateDistribution(muV,sigV);

% Stop simulations when we reach tmax = 10 minutes
Nevent = n0-1;          % array size only; run stops at tmax
tmin   = 2*60;          % window start [s]
tmax   = 10*60;         % window end [s]
xc     = mean(X0);      % centre of initial array [um]
hw     = 500;           % half-width of spatial window [um]

% Generate filename
muV_str  = strrep(num2str(muV, '%.3f'), '.', 'pt');
sigV_str  = strrep(num2str(sigV, '%.3f'), '.', 'pt');
outputfile = ['output_pooled_n_' num2str(n0) '_Nsim_' num2str(Nsim) '_muV_' muV_str '_sigV_' sigV_str '.mat']; 
save(outputfile,'Nsim','n0','dx0','X0','muV','sigV','vel_pdf','tmin','tmax','xc','hw')

% Simulations
obs_vel = cell(1,Nsim);
obs_len = cell(1,Nsim);
tic
for jj = 1:Nsim
    % Independent, reproducible random stream for each lane (Claude Opus 5.5 suggestion)
    % To get a fresh, independent set of simulations (for example, to check 
    % that a result isn’t a fluke of one seed), change 2026 to another number.
    s = RandStream('Threefry','Seed',2026);
    s.Substream = jj;
    RandStream.setGlobalStream(s);
    
    % Run simulation
    [x0,t0,tend,~,vel,len,way] = bactchainSimulation_pooled(n0,X0,vel_pdf,tmax,1,Nevent);

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

%% Analyze results (OPTIONAL)
nmax  = 8;                  % limit analysis to chains of length...
vel_all = [obs_vel{:}];
len_all = [obs_len{:}];

% Define variables
avg_vel = nan(1,nmax);
sem_vel = nan(1,nmax);

% Analyze results
for m = 1:nmax
    v = vel_all(len_all==m);
    avg_vel(m) = mean(v);
    sem_vel(m) = std(v)/sqrt(numel(v));
end

% Check 1: compare with experimental mean n = 2.2
meanLen     = mean(len_all);    
fprintf('Observed chains: %d, mean length: %.2f\n', numel(len_all), meanLen)

% Check 2: observed chains in each lane (1 x Nsim)
Nobs_lane = cellfun(@numel, obs_len);   
fprintf('Chains per lane: mean %.1f, min %d, max %d, total %d\n', ...
    mean(Nobs_lane), min(Nobs_lane), max(Nobs_lane), sum(Nobs_lane))

% Save results
save(outputfile,'nmax','avg_vel','sem_vel','meanLen','vel_all','len_all','obs_vel','obs_len','-append')