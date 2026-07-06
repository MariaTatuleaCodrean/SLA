%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
load('../Experiments/Codes/Output/chain_data.mat','expvel','explen')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define colors
col1 = [1 1 1];         % white
col2 = [178 30 51]/255; % burgundy
col3 = [22 134 118]/255;% dark turqouise
col4 = [30 178 158]/255;% turqouise
col5 = [0 0 0];         % black
col6 = [1 1 1]*.4;      % grey

edgcol = [1 1 1]*.5;

figure('Position',[530,535,248,137])

edges = 0.5:1:9.5;   % 
centers = (edges(1:end-1) + edges(2:end)) / 2;
barWidth = 0.8 * (centers(2) - centers(1));

bar(centers, histcounts(explen, edges), barWidth,'FaceColor',col6, 'FaceAlpha',0.5,'EdgeColor',edgcol);
xlabel('\textsf{chain length}, $n$','Interpreter','Latex')
ylabel('\textsf{count}','Interpreter','Latex')
set(gca, 'YTickLabelRotation', 0, 'FontSize', 9);
xlim([0.2 8.8])
xticks(1:1:8)
ylim([0 4500])
yticks([0 2000 4000])
box on
set(gca,'TickDir','out'); 
set(gca,'YAxisLocation', 'right');
legend('experiments')

% Save to file
print(gcf, 'Output/fig1e_inset.svg', '-dsvg');