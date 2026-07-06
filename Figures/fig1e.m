%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

% Load simulation & experimental results
load('../Experiments/Codes/Output/chain_data.mat','avg_v','sem_v')
outputfile = '../Simulations/Data/output_ABS_n_10000_Nsim_50_muV_9pt181_sigV_7pt255.mat'; 
load(outputfile,'dx0','Nsim','tplot','nmax','avg_vel','sigV','muV')

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
figure('Position', [750,438,313,238]);
hold on

% Plot marker for legend order
plot(-1,-1,'ko','MarkerFaceColor','k','DisplayName','experiments');

% Plot simulation predictions
x = 1:nmax;
tscale = dx0/muV;
y = avg_vel(:,:); %./avg_vel(1,:);

% plot each simulation as thin line
for kk=1:Nsim
    plot(x',y(:,kk),'-','Color',[mycols 0.1],'MarkerFaceColor',mycols,...
        'LineWidth',0.8,'HandleVisibility','off');
end

% plot mean as thick line
plot(x',mean(y,2),'o','Color',mycols,'MarkerFaceColor',mycols,...
    'DisplayName','simulations');

% Plot experimental measurements
x = 1:nmax;
y = avg_v; 
erry = sem_v; 
errorbar(x,y,erry,'ko','MarkerFaceColor','k','LineWidth',1.2,'HandleVisibility','off');

% x axis
xlim([.5 nmax+.5]);
xticks(1:nmax)
xlabel('\textsf{chain length,} $n$','Interpreter','Latex')

% y axis
% ylim([0.7 3.3])
% ylabel('\textsf{normalized speed,} $V_n/V_1$','Interpreter','Latex')
ylabel('\textsf{chain speed,} $V_n~\mathrm{[\mu m/s]}$','Interpreter','Latex')

% Manual northwest corner
legend('Location','southeast')

% Tidy up
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
box on
set(gca,'TickDir','out'); % The only other option is 'in'

% Title
% title(['\sigma_V = ' num2str(sigV) ' μm/s'])

% Save to file
print(gcf, 'Output/fig1e.svg', '-dsvg');