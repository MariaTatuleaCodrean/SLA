%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; addpath('functions') % load user defined functions

% Load experimental results
load('chain_data.mat','nvec','avg_v','sem_v','muV')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FITTING - Theoretical Model to Exp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess for sigV
sigV0 = muV;   % or choose a better guess if available

% Lower and upper bounds
lb = 0;
ub = Inf;

% Define model function for lsqcurvefit
modelFun = @(sigV, n) ABSFitFunction(n, muV, sigV);
% modelFun = @(sigV, n) ABSFitFunction_oneway(n, muV, sigV);
% modelFun = @(sigV, n) ABSGammaFitFunction(n, muV, sigV);

% Customize options for noisy model
options = optimoptions('lsqcurvefit', ...
    'FiniteDifferenceType','central',...
    'FiniteDifferenceStepSize',1e-1);

% Run nonlinear least squares fit
[ sigV_fit, resnorm, residual, exitflag, output ] = ...
    lsqcurvefit(modelFun, sigV0, nvec(:), avg_v(:), lb, ub,options);

% Reconstruct fitted curve
fitted_v = modelFun(sigV_fit, nvec);

% Save 
Nsim = 10;          % number of simulations
n0   = 1e4;         % number of bacteria 
dx0  = 1e3;         % regular spacing between bacteria (1mm)
sigV = sigV_fit;
save('modelfit_insertID.mat','Nsim','n0','dx0','muV','sigV')

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
plot(nvec, fitted_v, 'r-', 'LineWidth', 1, 'DisplayName', ['\sigma_V = ' num2str(sigV_fit, '%.2f')]);

% x axis
xlim([min(nvec)-.5 max(nvec)+.5]);
xticks(nvec)
xlabel('\textsf{chain length,} $n$','Interpreter','Latex')

% y axis
% ylim([0.7 3.3])
% ylabel('\textsf{normalized speed,} $V_n/V_1$','Interpreter','Latex')
ylabel('\textsf{chain speed,} $V_n$ [μm/s]','Interpreter','Latex')

% Manual northwest corner
legend('Location','northwest')

% Tidy up
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
box on
set(gca,'TickDir','out'); % The only other option is 'in'