%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% Load experimental results
csv_file = fullfile('../Experiments/Codes/Output/particle1_track_velocity.csv');  % Update filename if needed
data1 = readtable(csv_file);
T1 = data1.time_s;
V1 = data1.velocity_smooth_um_per_s;

csv_file = fullfile('../Experiments/Codes/Output/particle2_manual_velocity.csv');  % Update filename if needed
data2 = readtable(csv_file);
T2 = data2.time_s;
V2 = data2.velocity_smooth_um_per_s;

% Encounter time
tstar = T2(end);

% Chain
t = T1(T1>tstar) - tstar;
v = V1(T1>tstar);

% Particle 1
t1 = T1(T1<=tstar) - tstar;
v1 = V1(T1<=tstar);

% Particle 2
t2 = T2 - tstar;
v2 = V2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS - Experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift time so that merger is at t=0
merger_t = max(t1)/2 + min(t)/2; % time of merger
t  = t  - merger_t;
t1 = t1 - merger_t;
t2 = t2 - merger_t;

% Calculate average speeds before encounter
avg_v1 = mean(v1);
avg_v2 = mean(v2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
numcols = 3;
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
col5 = [0 0 0];         % black
mycols = colorgradient(numcols,(col4*0.6 + 0.4*col1),(col2*0.9 + 0.1*col1));

lnwdth = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT - Stacked Histogram Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [750,413,345,225]); %[750,438,345,200]);
hold on

% % Mean velocities of individual bacteria
plot([min(t1) max(t1)], [avg_v1 avg_v1],'--','Color',mycols(1,:),...
   'LineWidth',1.5,'HandleVisibility','off')
plot([min(t2) max(t2)], [avg_v2 avg_v2],'--','Color',mycols(3,:),...
   'LineWidth',1.5,'HandleVisibility','off')

% Velocity against time
plot(t1,v1,'-','Color',mycols(1,:),'LineWidth',lnwdth);
plot(t2,v2,'-','Color',mycols(3,:),'LineWidth',lnwdth); 
plot( t, v,'-','Color',mycols(2,:),'LineWidth',lnwdth,'HandleVisibility','off'); 

% chain velocity predicted
plot([0 7], [1 1]*(avg_v1+avg_v2)/2,'k--',...
   'LineWidth',1.5)

xlabel('\textsf{time,} $t ~\mathrm{[s]}$','Interpreter','Latex')
xlim([-40 7]);

ylabel('\textsf{speed,} $v ~\mathrm{[\mu m/s]}$','Interpreter','Latex')
ylim([0.2 4.2])
% yticks([0 0.05 0.10])
% yticklabels({'0','0.05','0.10'})

% Legend labels
leg_labels = {'rear bacterium', 'front bacterium', 'chain predicted'};

% Tidy up
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
box on
set(gca,'TickDir','out'); % The only other option is 'in'

% Save to file
print(gcf, 'Output/fig1d.svg', '-dsvg');