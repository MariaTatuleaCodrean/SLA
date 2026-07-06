% Preamble
clear; addpath('../Simulations/Codes/Functions') % load user defined functions
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

% Inputs
dx0 = 1;        % regular spacing between bacteria

% Generate initial speed distribution
muV  = 9;
sigV_vec = muV*[.3, 1, 2];
legend_labels = compose('%.1f', sigV_vec/muV);

% Assign control variables
Vc = muV*[.5, 1, 2];  % central bacterium velocity (for CDF plots)
rangeT = dx0/muV*linspace(0,10,2000); % range of times (for CDF plots)
rangeVc = muV*linspace(0.01,5,100);  % range of Vc (for median Te plots)

% Define output variables
all_cdfTe    = nan(length(Vc), length(rangeT),  length(sigV_vec));
all_pdfTe    = nan(length(Vc), length(rangeT),  length(sigV_vec));
all_Pinf     = nan(length(Vc), 1,               length(sigV_vec));
all_medianTe = nan(1,          length(rangeVc), length(sigV_vec));
all_finiteTe = nan(1,          length(rangeVc), length(sigV_vec));

%% Minimal model computations
tic

% For each sigV in sigV_vec
for kk=1:length(sigV_vec)
    % Compute encounter time statistics for each sigV 
    sigV = sigV_vec(kk);
    [~,mu,sigma] = generateDistribution(muV,sigV);
    
    % Generate output file name
    sigV_str  = strrep(num2str(sigV/muV, '%.3f'), '.', 'p');
    outputfile = ['../Simulations/Data/output_Trio_probability_theory_sigV_per_muV_' sigV_str '.mat'];

    [cdfTe, pdfTe, Pinf, medianTe, finiteTe] = bact3Simulation_Lognormal(mu,sigma,rangeVc,Vc,rangeT,dx0,outputfile);

    toc % time keeping

    all_cdfTe(:,:,kk)    = cdfTe; 
    all_pdfTe(:,:,kk)    = pdfTe; 
    all_Pinf(:,:,kk)     = Pinf; 
    all_medianTe(:,:,kk) = medianTe; 
    all_finiteTe(:,:,kk) = finiteTe; 
end

%% Post-processing
% load(outputfile,'cdfTe', 'Pinf', 'medianTe')
t = dx0/muV;       % characteristic time of simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
col5 = [0 0 0];         % black
edgcol = [1 1 1]*1;

% Define line styles
linestyle = {'-','--','-.',':'};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - Fig. 2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','centimeters','Position', [15,10,8.5,4]);
hold on

numcols = length(sigV_vec);
mycols = colorgradient(numcols,(col4*0.6 + 0.4*col1),(col2*0.9 + 0.1*col1));

for kk=1:length(sigV_vec)
    plot(rangeVc/muV,all_finiteTe(:,:,kk),linestyle{kk},'Color',mycols(kk,:),'LineWidth',1.75)
end

box on
set(gca,'TickDir','out'); % The only other option is 'in'

xlabel('\textsf{normalized speed,} $v_c/ \overline{V}$','Interpreter','Latex')
ylabel('$1 - P_{\infty}(v_c)$','Interpreter','Latex')

ylim([0.70 1.05])
yticks([0.75 1.0])
yticklabels({'3/4','1'})

ax = gca;
ax.FontSize = 9;

leg=legend(legend_labels,...
    'Location','southeast');
title(leg,'$\sigma_{V}/\overline{V}$','Interpreter','Latex')

print(gcf, 'Output/fig2b.svg', '-dsvg');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - Fig. 2c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','centimeters','Position', [15,10,8.5,5]);

% asymptotics
loglog([0.01 10],[1 1],'k--','LineWidth',1,'HandleVisibility','off')
hold on

numcols = length(sigV_vec);
mycols = colorgradient(numcols,(col4*0.6 + 0.4*col1),(col2*0.9 + 0.1*col1));

for kk=1:length(sigV_vec)
    loglog(rangeVc/muV,all_medianTe(:,:,kk)/t,'-','Color',mycols(kk,:),'LineWidth',1.75)
    hold on
end

% asymptotics
loglog(rangeVc/muV,       muV./rangeVc,'k--', 'LineWidth',1,'HandleVisibility','off')
ylim([0 10])

box on
set(gca,'TickDir','out'); % The only other option is 'in'

xlabel('\textsf{normalized speed,} $v_c/ \overline{V}$','Interpreter','Latex')
ylabel('\textsf{median,} $\widetilde{T}_e(v_c)/t_0$','Interpreter','Latex')

ax = gca;
ax.FontSize = 9;

leg=legend(legend_labels,...
    'Location','northeast');
title(leg,'$\sigma_{V}/\overline{V}$','Interpreter','Latex')

text(1.5*1e-2,3.5*1e-1,'$\widetilde{T}_e \to d_0/\overline{V}$','Interpreter','Latex')
text(1.65*1e-1,3.5*1e-1,'$\widetilde{T}_e \sim d_0/v_c$','Interpreter','Latex')


print(gcf, 'Output/fig2c.svg', '-dsvg');