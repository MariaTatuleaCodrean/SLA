%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
figpathname = 'Output/fig3h.svg';

mf = [1/3 2/3 1]; % multiplicative factor
tf = [3 2 1];     % tplot factor
legt = {'1/3','2/3','1'};  % legend labels

allt0   = {}; 
alltend = {}; 
alltsim = {}; 
alllen  = {};
allvel  = {};

for kk=1:length(mf)
    mf_str = strrep(num2str(mf(kk), '%.3f'), '.', 'p');
    outputfile = ['../Simulations/Data/output_ABS_n_10000_Nsim_10_sigL_per_muL_' mf_str 'Lognormal.mat'];

    load(outputfile,'dx0','t0','tend','tsim','len','vel','muV','Nsim')
    
    allt0{kk}   = t0;
    alltend{kk} = tend;
    alltsim{kk} = tsim;
    alllen{kk}  = max(1, round(len)); % round to nearest positive integer
    allvel{kk}  = vel;    
end

% Limit analysis to chains of length...
nmax = 8;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS - Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = dx0/muV;    % characteristic time of simulation
tplot = tf*t;  % time point for plotting

% Define variables
avg_vel = nan(length(mf),nmax,length(tplot),Nsim);
std_vel = nan(length(mf),nmax,length(tplot),Nsim);
sem_vel = nan(length(mf),nmax,length(tplot),Nsim);

% Determine evolution of chain length and speed
for kk=1:length(mf) % for each numerical experiment, ...
    for ii=1:Nsim % for each simulation, ...
        t0 = allt0{kk}(ii,:);
        tend = alltend{kk}(ii,:);
        tsim = alltsim{kk}(ii,:);
        len  = alllen{kk}(ii,:);
        vel  = allvel{kk}(ii,:);

        for jj=1:length(tplot) % at each time point, ...
            for m = 1:nmax % and each chain length 
                % Identify active chains of length m at time tkk
                I = find((t0<=tplot(jj)) & (tend>tplot(jj)) & (len==m));

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
col2 = [178 30 51]/255;     % burgundy
col3 = [22 134 118]/255;    % dark turqouise
col4 = [30 178 158]/255;    % turqouise
col5 = [0 0 0];             % black
col6 = [107 92 165]/255;    % purple
col7 = [242, 201 76]/255;   % yellow
mycols = colorgradient(length(mf),(col7*0.9 + 0.1*col1),(col6*0.8 + 0.2*col1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure('Units', 'centimeters','Position', [15,10,5.6,5.6]);
hold on;

% Plot simulation predictions
x = 1:nmax;
for jj=1:length(mf)
    mypastelcols = colorgradient(length(tplot),mycols(jj,:),col1);

    for ii = [3 1]
        if ii==1
            hvis = 'on';
        else
            hvis = 'off';
        end
        y = avg_vel(jj,:,ii,:)/muV;
        plot(x,mean(y,4),'o-','Color',mycols(jj,:),...
            'MarkerFaceColor',mypastelcols(ii,:),...
            'HandleVisibility',hvis,...
            'DisplayName',legt{jj});
    end
end

% x axis
xlim([.5 8.5]);
xticks(1:nmax)
xlabel('\textsf{chain length,} $\hat{n}(L)$','Interpreter','Latex')

% y axis
ylim([0.5 2.5])
ylabel('\textsf{mean speed,} $V_n/\overline{V}$','Interpreter','Latex')

% Legend
leg = legend('Location','northwest','FontSize',7);
title(leg,'$\sigma_L/ \overline{L}$','Interpreter','Latex')

% Tidy up
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
box on
set(gca,'TickDir','out'); 

% Save to file
print(gcf, figpathname, '-dsvg');