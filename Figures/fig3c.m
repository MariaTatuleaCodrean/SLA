%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

load('../Simulations/Data/output_ABS_n_10000_Nsim_50_muV_9pt000_sigV_6pt000.mat',...
    'dx0','vel_pdf','Nsim','t0','tend','tsim','vel','len','muV') 
figpathname = 'Output/fig3c.svg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = dx0/muV;                        % characteristic time of simulation
tplot = linspace(0,3,100)*t;        % time range for plotting
% tplot = setdiff(tplot,tsim(:));   % remove times of merging events

length_val = (1:1:7)';   % analyze chains with these lengths

% Define plot variables
length_count = nan(length(length_val),length(tplot),Nsim);
average_vel  = nan(length(length_val),length(tplot),Nsim);

% Determine evolution of chain length and speed
for jj=1:Nsim % for each simulation, ...
    for kk=1:length(tplot) % at each time point, ...
        tkk = tplot(kk);
        for ii = 1:length(length_val) % and each chain length
            m = length_val(ii);

            % Identify active chains of length m at time tkk
            I = find((t0(jj,:)<=tkk) & (tend(jj,:)>tkk) & (len(jj,:)==m));
            length_count(ii,kk,jj) = size(I,2);
            average_vel(ii,kk,jj)  = mean(vel(jj,I));
        end
    end
end

% X values
X = tplot/t; % common time axis

% Y values | Compute statistics across the Nsim simulations
Y_mean = mean(average_vel,3); 
Y_std  = std(average_vel,0,3); % standard deviation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
numcols = size(Y_mean,1);
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
mycols = colorgradient(numcols,(col2*0.8 + 0.2*col1),(col4*0.8 + 0.2*col1));

edgcol = [1 1 1]*1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - overlay plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure('Units', 'centimeters','Position', [15,10,8.26,5.3]); 
hold on;

% Plot chains of length 8, 8--7, 8--6, ..., until 8--1
for kk = 1:size(Y_mean,1)
    mycol = mycols(numcols+1-kk,:);

    % Every simulation as thin line
    for jj=1:size(average_vel,3)
        plot(X,average_vel(kk,:,jj), 'Color', [mycol 0.1], 'LineWidth', 0.2,...
            'HandleVisibility','off');
    end

    % Mean as thick line
    plot(X,Y_mean(kk,:), 'Color', mycol, 'LineWidth', 1.2,'DisplayName',['n = ' num2str(length_val(kk))]);
end

xlim([0 max(X)]);
ylim([0 50]);
box on
xticks([0 1 2 3])
set(gca,'TickDir','out'); % The only other option is 'in'

xlabel('\textsf{normalized time,} $t/t_0$','Interpreter','Latex')
ylabel('\textsf{mean speed,} $V_n~\mathrm{[\mu m/s]}$','Interpreter','Latex')

ax = gca; 
ax.FontSize = 9; 

[h,icons]=legend('Position',[0.662728721920143,0.379882546556881,0.236760122289539,0.524999985396862]);
set(h,'FontSize',7);
legend('boxoff') %This is optional
h.ItemTokenSize = [10,15];

% Save to file
print(gcf, figpathname, '-dsvg');