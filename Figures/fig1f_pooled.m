%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

% Load experimental results
load('chain_data.mat','expvel','explen')

% Load simulation results
outputfile = 'output_pooled_n_20_Nsim_1000_muV_9pt181_sigV_10pt323.mat'; % identify simulation
load(outputfile,'len_all','vel_all','Nsim','vel_pdf')

% Only plot selected n values
selection = [1 3 5];
ctmylim = [0.15 0.1 0.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
numcols = length(selection);
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
col5 = [0 0 0];         % black
mycols = colorgradient(numcols,(col4*0.6 + 0.4*col1),(col2*0.9 + 0.1*col1));
mycols2 = colorgradient(numcols,(col4*0.6 + 0.4*col5),(col2*0.9 + 0.1*col5));

edgcol = [1 1 1]*0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - Stacked Histogram Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [750,324,313,332]); % 3 rows & 1 column %329

% Loop through subplots
for kk = 1:length(selection)
    % Focus on m-length chains
    m = selection(kk);

    subplot(3,1,kk)
    hold on

    % Histogram plot - experimental data
    histogram(expvel(explen == m),'Normalization','pdf','BinWidth',1, ...
        'FaceColor',mycols(kk,:),'EdgeColor',edgcol);
 
    % Fit velocity distribution from simulations
    all_vel = vel_all(len_all==m);
    sim_pdf = fitdist(all_vel(:), 'Kernel', 'Support', 'positive');   % fit kernel density to x

    % Plot
    x_values = linspace(0,40);
    y_values = pdf(sim_pdf,x_values);
    plot(x_values,y_values,'Color',[mycols2(kk,:) 1],'LineWidth',1.5,'HandleVisibility','on')

    % Tidy up
    xlim([0 40]);
    % ylim([0 0.3])
    yticks(0:0.1:0.3)
    yticklabels({'0','0.1','0.2','0.3'})  

    box on
    set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
    set(gca,'TickDir','out'); % The only other option is 'in'

    if kk==1
        % Initial distribution for first plot
        x_values = linspace(0,40);
        y_values = pdf(vel_pdf,x_values);
        plot(x_values,y_values,'k-','LineWidth',1.5,'HandleVisibility','on')

        legend('\textsf{experiments}','\textsf{theory (pooled)}','\textsf{theory} ($t=0$)','Interpreter','Latex')
    else
        legend('\textsf{experiments}','\textsf{theory (pooled)}','Interpreter','Latex')
    end

    lgd = legend;      % get current legend handle
    lgd.FontSize = 8; % change font size

    if kk==length(selection)
        % x-label for last subplot
        xlabel('\textsf{speed,} $v~\mathrm{[\mu m/s]}$','Interpreter','Latex')
    end
    switch kk
        case 1
            % y-label for middle subplot
            ylabel('\textsf{PDF,} $f_{V}(v)$','Interpreter','Latex')
        case 2
            % y-label for middle subplot
            ylabel('\textsf{PDF,} $f_{V}(v)$','Interpreter','Latex')
        case 3
            % y-label for middle subplot
            ylabel('\textsf{PDF,} $f_{V}(v)$','Interpreter','Latex')
    end
    % ylabel('\textsf{PDF}','Interpreter','Latex')

    % xlabel('\textsf{speed,} $v~\mathrm{[\mu m/s]}$','Interpreter','Latex')
    % ylabel('\textsf{PDF,} $f_{V}(v)$','Interpreter','Latex')

    % Add label for chain length
    % txt = ['n = ' num2str(selection(kk)) ' (N = ' num2str(numel(expvel)) ')'];
    % text(20, 0.115,txt,'FontWeight','bold')
    % title(txt)    
end

% Save to file
print(gcf, 'figures\fig1f_revised_pooled.svg', '-dsvg');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - Chain length distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Longest chain length shown (longer chains still count in the normalization)
nmax = 8;
 
% Fraction of observed chains with each length n = 1..nmax
% (normalized by ALL chains, so fractions are comparable between data sets)
edges   = 0.5:1:nmax+0.5;
centers = 1:nmax;
Nexp = numel(explen);
Nsim_chains = numel(len_all);
f_exp = histcounts(explen,  edges) / Nexp;
f_sim = histcounts(len_all, edges) / Nsim_chains;
 
% Binomial standard error of the experimental fractions
err_exp = sqrt(f_exp.*(1-f_exp)/Nexp);
 
% Mean chain length (all chains, not just n <= nmax)
meanLen_exp = mean(explen);
meanLen_sim = mean(len_all);
fprintf('Mean chain length: experiments %.2f (N = %d), simulations %.2f (N = %d)\n', ...
    meanLen_exp, Nexp, meanLen_sim, Nsim_chains)
 
% Make new figure
figure('Position', [750,438,313,238]);
hold on
 
% Side-by-side bars
barWidth = 0.4;
bar(centers - barWidth/2, f_exp, barWidth, 'FaceColor',col2, 'FaceAlpha',0.6, 'EdgeColor',edgcol, ...
    'DisplayName', sprintf('experiments, \\langle n\\rangle = %.2f', meanLen_exp));
bar(centers + barWidth/2, f_sim, barWidth, 'FaceColor',col4, 'FaceAlpha',0.6, 'EdgeColor',edgcol, ...
    'DisplayName', sprintf('theory (pooled), \\langle n\\rangle = %.2f', meanLen_sim));
 
% Error bars on the experimental bars
errorbar(centers - barWidth/2, f_exp, err_exp, 'k', 'LineStyle','none', 'HandleVisibility','off');
 
% x axis
xlim([0.5 nmax+0.5]);
xticks(1:nmax)
xlabel('chain length, {\itn}')
 
% y axis
ylabel('fraction of chains')
 
% Legend
lgd = legend('Location','northeast');
lgd.FontSize = 8;
legend boxoff
 
% Tidy up
box on
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
set(gca,'TickDir','out'); % The only other option is 'in'