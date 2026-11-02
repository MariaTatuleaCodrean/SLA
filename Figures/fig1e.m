%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

% Load simulation & experimental results
load('../Experiments/Codes/Output/chain_data.mat','vel_pdf_V1','explen')
outputfile = '../Simulations/Data/modelfit_29October2026.mat'; 
load(outputfile,'avg_v', 'sem_v', 'nvec', 'fitted_v', 'fitted_sem_v','sigV')

disp(['best fit sigV = ' num2str(sigV)])
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM ASSEMBLY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges = 0.5:1:8.5;
Nvec = histcounts(explen, edges); % count number of chains of each length
random_avgV = nan(1,numel(Nvec));  % rows: sample, columns: chain length
random_semV = nan(1,numel(Nvec));  % rows: sample, columns: chain length
for kk = 1:numel(Nvec)
    % Generate randomly assembled chains
    random_bactrV = random(vel_pdf_V1,[kk Nvec(kk)]);
    random_chainV = mean(random_bactrV, 1);

    % Mean chain speed for this length in this sample
    random_avgV(kk) = mean(random_chainV);
    random_semV(kk) = std(random_chainV)/sqrt(numel(random_chainV));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
col5 = [0 0 0];         % black

numcols = length(outputfile);
% mycols = colorgradient(numcols,(col4*0.6 + 0.4*col1),(col2*0.9 + 0.1*col1));
mycols = col4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [750,438,313,268]);
hold on

% Plot marker for legend order
plot(-1,-1,'ko','MarkerFaceColor','k','DisplayName','experiments');

% Plot experimental measurements
errorbar(nvec,avg_v,sem_v,'ko','MarkerFaceColor','k','LineWidth',1.2,'HandleVisibility','off');

% Plot simulation predictions
errorbar(nvec, fitted_v, fitted_sem_v, 'o', 'MarkerFaceColor', col4, 'Color', col4, ...
    'LineWidth', 1, 'DisplayName', 'simulations (pooled, 2-10 min)');

% Plot marker for legend order
plot(-1,-1,'o','Color',col2,'MarkerFaceColor',col2,'DisplayName','theory – random assembly');

% Plot random predictions: each sample as thin line
errorbar(nvec,random_avgV,random_semV,'o','Color',col2,'MarkerFaceColor',col2,'HandleVisibility','off'); % mean across samples

% x axis
xlim([.5 max(nvec)+.5]);
xticks(1:max(nvec))
xlabel('\textsf{chain length,} $n$','Interpreter','Latex')

% y axis
ylim([2 21])
% ylabel('\textsf{normalized speed,} $V_n/V_1$','Interpreter','Latex')
ylabel('\textsf{chain speed,} $V_n~\mathrm{[\mu m/s]}$','Interpreter','Latex')

% Manual northwest corner
% legend('Location','southeast')

% Tidy up
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
box on
set(gca,'TickDir','out'); % The only other option is 'in'

% Title
% title(['\sigma_V = ' num2str(sigV) ' μm/s'])

% Save to file
print(gcf, 'Output/fig1e.svg', '-dsvg');