%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3e WITH INSET
% Main plot: mean chain length vs normalized time
% Inset: mean speed vs chain length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

figpathname = 'Output/fig3e.svg';

Pright = [0.5 0.8 1.0];         % right-moving probability
legt  = {'0.5','0.8','1.0'};    % legend labels

tmax = Inf;                     % limit analysis to available times
nmax_inset = 8;                 % maximum chain length in inset
tf_inset = [3 2 1];             % tplot factors used for inset

allt0   = {}; 
alltend = {}; 
alltsim = {}; 
alllen  = {};
allvel  = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:length(Pright)
    P_str = strrep(num2str(Pright(kk), '%.2f'), '.', 'p');
    outputfile = ['../Simulations/Data/output_ABS_n_10000_Nsim_10_Pright_' P_str '.mat'];

    load(outputfile,'dx0','t0','tend','tsim','len','vel','muV','Nsim')
    
    allt0{kk}   = t0;
    alltend{kk} = tend;
    alltsim{kk} = tsim;
    alllen{kk}  = len;
    allvel{kk}  = vel;
    
    % Limit analysis to shortest simulation time
    tmax = min(tmax, min(tsim(:,end)));
end

% Characteristic time of simulation
tchar = dx0/muV;

% Limit main-plot analysis to chain lengths observed in simulations
nmax_main = max([alllen{:}],[],'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS MAIN PLOT: mean chain length vs time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tplot_main = linspace(0,tmax,100);

% Define variables
length_count = zeros(length(Pright),length(tplot_main),Nsim);
agent_count  = zeros(length(Pright),length(tplot_main),Nsim);
avg_len      = zeros(length(Pright),length(tplot_main),Nsim);

% Determine evolution of chain length and speed
for kk = 1:length(Pright) % for each numerical experiment, ...
    for ii = 1:Nsim % for each simulation, ...
        t0   = allt0{kk}(ii,:);
        tend = alltend{kk}(ii,:);
        len  = alllen{kk}(ii,:);

        for jj = 1:length(tplot_main)
            for m = 1:nmax_main
                % Identify active chains of length m at time tkk
                I = find((t0 <= tplot_main(jj)) & (tend > tplot_main(jj)) & (len == m));
                length_count(kk,jj,ii) = length_count(kk,jj,ii) + m*size(I,2);
                agent_count(kk,jj,ii)  = agent_count(kk,jj,ii)  +   size(I,2);
            end
        end

        avg_len(kk,:,ii) = length_count(kk,:,ii)./agent_count(kk,:,ii);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS INSET: mean speed vs chain length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tplot_inset = tf_inset*tchar; % time points for plotting

avg_vel = nan(length(Pright),nmax_inset,length(tplot_inset),Nsim);

% Determine evolution of chain length and speed
for kk = 1:length(Pright) % for each numerical experiment, ...
    for ii = 1:Nsim % for each simulation, ...
        t0   = allt0{kk}(ii,:);
        tend = alltend{kk}(ii,:);
        len  = alllen{kk}(ii,:);
        vel  = allvel{kk}(ii,:);

        for jj = 1:length(tplot_inset) % at each time point, ...
            for m = 1:nmax_inset % and each chain length 
                % Identify active chains of length m at time tkk
                I = find((t0 <= tplot_inset(jj)) & (tend > tplot_inset(jj)) & (len == m));

                % Calculate average velocity
                avg_vel(kk,m,jj,ii) = mean(vel(I));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
col1 = [1 1 1];             % white
col6 = [107 92 165]/255;    % purple
col7 = [242, 201 76]/255;   % yellow
mycols = colorgradient(length(Pright),(col7*0.9 + 0.1*col1),(col6*0.8 + 0.2*col1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT MAIN FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Units', 'centimeters','Position', [15,10,8.27,5.5]);

ax_main = axes(fig);
hold(ax_main,'on')

x_main = tplot_main/tchar;
for kk = 1:length(Pright)
    y = avg_len(kk,:,:);
    plot(ax_main,x_main,mean(y,3),'-','Color',mycols(kk,:),...
        'LineWidth',1.5,...
        'MarkerFaceColor',mycols(kk,:),...
        'DisplayName',legt{kk});
end

xlim(ax_main,[0 3]);
xlabel(ax_main,'\textsf{normalized time,} $t/t_0$','Interpreter','Latex')

ylim(ax_main,[0.7 6.7])
ylabel(ax_main,'\textsf{mean chain length,} $\langle n \rangle$','Interpreter','Latex')

leg = legend(ax_main,'Location','northeast','FontSize',7);
title(leg,'$P(\textrm{right})$','Interpreter','Latex')

set(ax_main, 'YTickLabelRotation', 0, 'FontSize', 9);
box(ax_main,'on')
set(ax_main,'TickDir','out')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT INSET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position is [left bottom width height] in normalized figure units.
% Adjust these four numbers to move/resize the inset.
inset_pos = [0.15,0.45,0.31,0.44];
ax_inset = axes(fig,'Position',inset_pos);
hold(ax_inset,'on')

x_inset = 1:nmax_inset;
for kk = 1:length(Pright)
    mypastelcols = colorgradient(length(tplot_inset),mycols(kk,:),col1);

    % Plot only jj = 2, i.e. tplot_inset(2) = 2*tchar.
    jj = 2;
    y = avg_vel(kk,:,jj,:)/muV;
    plot(ax_inset,x_inset,mean(y,4),'o-','Color',mycols(kk,:),...
        'MarkerSize',5,...
        'LineWidth',0.5,...
        'MarkerFaceColor',mypastelcols(jj,:),...
        'HandleVisibility','off');
end

xlim(ax_inset,[0.5 8.5]);
xticks(ax_inset,1:nmax_inset)
ylim(ax_inset,[0.5 2])
yticks([1 2])

xlabel(ax_inset,'$n$','Interpreter','Latex','FontSize',9)
ylabel(ax_inset,'$V_n/\overline{V}$','Interpreter','Latex','FontSize',9)

set(ax_inset, 'YAxisLocation', 'right', 'XTickLabelRotation', 0) %'XAxisLocation', 'top',
set(ax_inset,'FontSize',9,'TickDir','out','TickLength',[0.02, 0.01])
box(ax_inset,'on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE TO FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(fig, figpathname, '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generated by ChatGPT (Model GPT-5.5 Thinking) & verified by user. 
% AI prompt: 
% Hello. I'd like your help with making some figures in MATLAB. 
% I will give you two separate MATLAB scripts, and I'd like you to combine 
% them so that one plot appears as an inset in the main figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%