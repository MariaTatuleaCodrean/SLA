%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

load('../Simulations/Data/output_ABS_n_10000_Nsim_50_muV_9pt000_sigV_6pt000.mat',...
    'n0','dx0','vel_pdf','Nsim','t0','tend','tsim','vel','len','muV') 
simk = 1; % choose specific simulation for this plot
figpathname = 'Output/fig3b.svg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = dx0/muV;                         % characteristic time of simulation
tplot = linspace(0,3,1000)*t;        % time range for plotting
tplot = setdiff(tplot,tsim(simk,:)); % remove times of merging events
tend(isnan(tend)) = Inf;             % include chains active at the end

length_val = (1:1:8)';   % analyze chains with these lengths

% Define plot variables
length_count = nan(length(length_val),length(tplot));
average_vel  = nan(length(length_val),length(tplot));

% Determine evolution of chain length and speed
length_val = (1:1:8)';
for kk=1:size(tplot,2)
    tkk = tplot(kk);
    for ii = 1:length(length_val)
        m = length_val(ii);

        if ii < length(length_val)
            % Identify active chains of length m at time tkk
            I = find((t0(simk,:)<tkk) & (tend(simk,:)>tkk) & (len(simk,:) ==m));
        else
            % Identify active chains of length 8+ at time tkk
            I = find((t0(simk,:)<tkk) & (tend(simk,:)>tkk) & (len(simk,:) >=m));
        end
        length_count(m,kk) = size(I,2); 
        average_vel(m,kk)  = mean(vel(I));
    end
end

% Flip arrays so chain length goes from 8 to 1
length_val   = flip(length_val,1);
length_count = flip(length_count,1);
average_vel  = flip(average_vel,1);

% X values
X = tplot/t; 

% Y values
Y = cumsum(length_count,1); % cumulative density from chain length m to 1
Y0 = zeros(1, length(X));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
numcols = size(Y,1);
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
mycols = colorgradient(numcols,(col2*0.7 + 0.3*col1),(col4*0.2 + 0.8*col1));

edgcol = [1 1 1]*1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - Stacked Band Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure('Units', 'centimeters','Position', [15,10,8.5,5.4]);
hold on;

% Plot chains of length 8, 8--7, 8--6, ..., until 8--1
for kk = size(Y,1):-1:1 
    % Fill between Y and Y0
    fill([X fliplr(X)], [Y(kk,:) fliplr(Y0)], mycols(kk,:), 'EdgeColor', edgcol, 'LineWidth', 0.25);
end

% Plot tstar, Nstar point
temp = Y(end,:);
I = find(temp<0.45*n0,1);
plot(X(I),temp(I),'ko','MarkerSize',4,'MarkerFaceColor','k')

xlim([0 max(tplot/t)]);
ylim([0 1e4]);
box on
xticks([0 1 2 3])
yticks([0 .25 .5 .75 1]*1e4)
ax = gca;
ax.YAxis.Exponent = 4;   % for x-axis
set(gca,'TickDir','out'); % The only other option is 'in'

xlabel('\textsf{normalized time,} $t/t_0$','Interpreter','Latex')
ylabel('\textsf{chain count,} $N_n(t)$','Interpreter','Latex')

ax = gca; 
ax.FontSize = 9; 

[h,icons]=legend({' n = 1',' n = 2',' n = 3',' n = 4',' n = 5',' n = 6',' n = 7',' n \geq 8'},...
    'Position',[0.634708976824537,0.327413993991831,0.24299065169887,0.557499984502793]);
% set(h,'NumColumns',2);
set(h,'FontSize',7);
legend('boxoff') %This is optional

boxes=findobj(icons,'type','patch');
for kk=1:8
    boxx=boxes(kk);
    boxx.XData=[0.4;0.4;0.6;0.6];
    boxx.LineWidth=1;
end

% Save to file
print(gcf, figpathname, '-dsvg');