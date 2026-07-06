%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA - Experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
csv_file = fullfile('../Data/chain_data.csv');  % Update filename if needed
data = readtable(csv_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS - Individual Bacteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select individual bacteria
data1 = data(data.chain_length == 1, :);

% Calculate average velocity of single bacteria
V1 = mean(data1.velocity);

% Fit velocity distributions
vel_pdf_V1 = fitdist(data1.velocity, 'Lognormal');   % fit lognormal density to x
mu_V1 = vel_pdf_V1.mu;
sig_V1 = vel_pdf_V1.sigma;

% Smooth pdf
x_values = linspace(0,40);
y = pdf(vel_pdf_V1,x_values);

% Create figure
figure;
hold on;

% Plot histogram for experimental data & distribution fit
histogram(data1.velocity, 'Normalization', 'pdf', 'DisplayName', 'Experimental sample');
plot(x_values,y,'r-','LineWidth',2,'DisplayName','Lognormal fit')

legend;
hold off;

box on
xlim([0 40])
ylabel('probability density')
xlabel('velocity [μm s^{-1}]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS - All Chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All velocities & all lengths
[expvel, explen] = deal(data.velocity, data.chain_length); 
[expvel, explen] = deal(expvel(:), explen(:));

% Limit analysis to chains of length n=8
nvec = 1:8;

% Define variables
avg_v = nan(size(nvec));
std_v = nan(size(nvec));
sem_v = nan(size(nvec));

% Loop through chain length
for kk = 1:length(nvec)
    % Select chains of length n
    data_kk = data(data.chain_length == nvec(kk), :);
    
    % Calculate average velocity
    avg_v(kk) = mean(data_kk.velocity);

    % Calculate standard deviation
    w = 0; %  the standard deviation is normalized by N-1
    std_v(kk) = std(data_kk.velocity,w); 

    % Calculate standard error of the mean
    % Recall that s.e.m. = std(x)/sqrt(length(x));
    sem_v(kk) = std_v(kk)/sqrt(length(data_kk.velocity));
end

% Speed (normalized by V1)
norm_v = avg_v/V1;

% Calculate number of tracked agents and total number of bacteria
agentsN = length(data.chain_length);
bacterN = sum(data.chain_length);

% Calculate active fraction in experiments = active agents/total # bacteria
alpha_exp = agentsN/bacterN;

% Calculate average number of bacteria per chain
avg_n = bacterN/agentsN;

% Calculate true mean speed statistic
muV = sum((data.velocity).*(data.chain_length))/bacterN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Output/chain_data.mat','data','V1','vel_pdf_V1','mu_V1','sig_V1',...
    'nvec','avg_v','std_v','sem_v','norm_v','muV',...
    'agentsN','bacterN','alpha_exp','avg_n','expvel','explen')