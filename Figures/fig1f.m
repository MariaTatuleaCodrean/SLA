%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

% Load experimental results
load('../Experiments/Codes/Output/chain_data.mat','expvel','explen')

% Load simulation results
outputfile = '../Simulations/Data/output_ABS_n_10000_Nsim_50_muV_9pt181_sigV_7pt255.mat'; 
load(outputfile,'sim_vel','sim_len','Nsim','vel_pdf')

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
% figure('Position', [750,291,313,385]); % 3 rows & 1 column (taller)
figure('Position', [750,324,313,352]); % 3 rows & 1 column %329
% figure('Position', [416,137,678,557]); % 4 rows & 2 columns

% Loop through subplots
for kk = 1:length(selection)
    % Focus on m-length chains
    m = selection(kk);

    subplot(3,1,kk)
    hold on

    % Histogram plot
    histogram(expvel(explen == m),'Normalization','pdf','BinWidth',1, ...
        'FaceColor',mycols(kk,:),'EdgeColor',edgcol);

    % % Add line for mean
    % plot([mean(expvel(explen == m)) mean(expvel(explen == m))], [0 0.15],'k','Color',mycols2(kk,:),...
    %     'LineWidth',1.5,'HandleVisibility','off') 
    
    for jj=1:Nsim
        % Extract current simulation data
        vel = sim_vel{jj};
        len = sim_len{jj};
        all_vel = vel(round(len)==m);

        if ~isempty(all_vel)
            % Fit velocity distribution from simulations
            sim_pdf = fitdist(all_vel(:), 'Kernel', 'Support', 'positive');   % fit kernel density to x

            % Plot
            x_values = linspace(0,40);
            y_values = pdf(sim_pdf,x_values);

            plot(x_values,y_values,'Color',[mycols2(kk,:) 0.1],'LineWidth',1,'HandleVisibility','off')
        end
    end

    % Draw hidden line for legend
    plot(-x_values,y_values,'Color',[mycols2(kk,:) 1],'LineWidth',1.5,'HandleVisibility','on')


    % Tidy up
    xlim([0 35]);
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

        legend('\textsf{experiments}','\textsf{theory} ($t=t_*$)','\textsf{theory} ($t=0$)','Interpreter','Latex')
    else
        legend('\textsf{experiments}','\textsf{theory} ($t=t_*$)','Interpreter','Latex')
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

    % Add label for chain length
    % txt = ['n = ' num2str(selection(kk)) ' (N = ' num2str(numel(expvel)) ')'];
    % text(20, 0.115,txt,'FontWeight','bold')
    % title(txt)    
end

% Save to file
print(gcf, 'Output/fig1f.svg', '-dsvg');