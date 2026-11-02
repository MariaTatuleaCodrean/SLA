%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; addpath('functions') % load user defined functions

% Load experimental results
load('../../Experiments/Codes/Output/chain_data.mat','nvec','avg_v','sem_v','muV')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FITTING - Theoretical Model to Experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
n0      = 20;          % number of bacteria per lane 
Nsim    = 1000;        % number of simulations
tmin    = 2*60;        % window start [s]
tmax    = 10*60;       % window end [s]
hw      = 500;         % half-width of spatial window [um]
myseed  = 2026;        % default seed for reproducible random initial conditions

% Scan range
nGridPoints = 31;
sigVvec = linspace(5.3436, 19.0843, nGridPoints);

% Initialize output variables
sim_v   = nan(nGridPoints, numel(nvec));
sim_sem = nan(nGridPoints, numel(nvec));
sim_N   = nan(nGridPoints, 1);
chi2    = nan(1, nGridPoints);

% Simulate system across chosen sigV range
outputfile = 'modelfit_insertID.mat';
for kk = 1:nGridPoints
    [sim_v(kk,:), sim_sem(kk,:), sim_N(kk)] = ABSFitFunction_pooled(nvec,muV,sigVvec(kk),n0,Nsim,tmin,tmax,hw,myseed);

    % Welch's t-test (adapted): residuals weighted by combined uncertainty
    t = (sim_v(kk,:) - avg_v) ./ sqrt(sem_v.^2 + sim_sem(kk,:).^2);
    chi2(kk) = sum(t.^2);
   
    % Print diagnostics and save progress after every grid point
    fprintf('sigV = %.2f um/s: chi2 = %.2f\n', sigVvec(kk), chi2(kk))
    save(outputfile,'sigVvec','sim_v','sim_sem','chi2','nvec','avg_v','sem_v','muV','n0','Nsim','sim_N')
end

% Diagnostic: Check total number of detected chains
disp(sim_N)

% Find lowest error
[~, ib] = min(chi2);
if ib == 1 || ib == nGridPoints
    warning('Minimum at the edge of the grid: widen sigGrid and rerun.')
end

% Fit parabola around minimum turning point
nb = max(1,ib-2):min(nGridPoints,ib+2);     % up to 5 points around minimum
p  = polyfit(sigVvec(nb), chi2(nb), 2);     % chi2 ~ a*(s-s0)^2 + c

% Best fit sigV is the minimum of the fitted parabola
sigV_fit = -p(2)/(2*p(1));
fprintf('sigV = %.2f um/s\n', sigV_fit)

% Diagnostic plot: Check that error looks like a parabola
figure; hold on
plot(sigVvec, chi2,'bo-')
plot(sigV_fit, p(1)*sigV_fit^2 + p(2)*sigV_fit + p(3),'ro','MarkerFaceColor','r')
xlabel('\sigma_V [\mu m/s]'); ylabel('\chi^2')

% Reconstruct fitted curve
[fitted_v, fitted_sem_v, fitted_sim_N] = ABSFitFunction_pooled(nvec,muV,sigV_fit,n0,Nsim,tmin,tmax,hw);

% Save 
sigV = sigV_fit;
save(outputfile,'sigV','fitted_v','fitted_sem_v','tmin','tmax','fitted_sim_N','-append')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
col5 = [0 0 0];         % black

edgcol = [1 1 1]*1;

% Make new figure
figure('Position', [750,438,313,238]);
hold on

% Plot experimental measurements
errorbar(nvec, avg_v, sem_v, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'DisplayName', 'experiments');

% Plot fitted curve
errorbar(nvec, fitted_v, fitted_sem_v, 'o', 'MarkerFaceColor', col4, 'Color', col4, ...
    'LineWidth', 1, 'DisplayName', 'simulations (pooled, 2-10 min)');

% x axis
xlim([min(nvec)-.5 max(nvec)+.5]);
xticks(nvec)
xlabel('\textsf{chain length,} $n$','Interpreter','Latex')

% y axis
ylabel('\textsf{chain speed,} $V_n~\mathrm{[\mu m/s]}$','Interpreter','Latex')

% Manual northwest corner
legend('Location','northwest')

% Tidy up
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
box on
set(gca,'TickDir','out'); % The only other option is 'in'