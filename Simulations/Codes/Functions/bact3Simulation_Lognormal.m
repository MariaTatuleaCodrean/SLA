function [cdfTe, pdfTe, Pinf, medianTe, finiteTe] = bact3Simulation_Lognormal(mu,sigma,rangeVc,Vc,rangeT,dx0,outputfile)
%BACT3SIMULATION Simulation of bacterial trio on single lane swimming to 
% predict distribution of encounter time and probability of finite-time
% encounter.
%
% Input arg     Type            Description
%   mu          double          mean of Lognormal distribution for speed
%   sigma       double          SD of Lognormal distribution for speed
%   rangeVc     double array    central bacterium velocity (for median Te)
%   Vc          double array    central bacterium velocity (for CDF plot)
%   rangeT      double array    time (for CDF plot)
%   dx0         double          initial spacing between bacteria (regular)
% 
% Optional arguments
%   outputfile  str         name of outputfile   
%
% Output arg
%   cdfTe       double array    cdfTe(j,k) = P(Te < rangeT(k) | Vc(j))
%   pdfTe       double array    probability density function | Vc(j)  
%   Pinf        double array    infinite-time encounter prob. | Vc(j)
%   medianTe    double array    median encounter time | rangeVc(j)
%   finiteTe    double array    finite-time encounter prob. | rangeVc(j)

% Define output variables
cdfTe = nan(length(Vc),length(rangeT));
pdfTe = nan(length(Vc),length(rangeT));
Pinf  = nan(length(Vc),1);
finiteTe = nan(1,length(rangeVc));
medianTe = nan(1,length(rangeVc));

% Compute mean speed
muV = exp(mu+sigma^2/2);

% For each vc in Vc
for jj=1:length(Vc)
    % Velocity of central bacterium
    vc = Vc(jj);

    % Evaluate CDF of encounter time 
    cdfTe(jj,:) = computeCDFTe_Lognormal(mu,sigma,vc,dx0,rangeT);

    % Evalute PDF of encounter time
    pdfTe(jj,:) = gradient(cdfTe(jj,:), rangeT);
    
    % Probability of infinite encounter time
    Fv_vc = logncdf(vc,mu,sigma);
    Pinf(jj) = (1-Fv_vc^2)/4; % always < 1/4, so median time is well-defined

    % Normalize so that integral over rangeT equals 1-Pinf
    pdfTe(jj,:) = (1-Pinf(jj)) * pdfTe(jj,:) ./ trapz(rangeT, pdfTe(jj,:));
end

% For each vc in rangeVc
for kk=1:length(rangeVc)
    % Velocity of central bacterium
    vc = rangeVc(kk);

    % Probability of infinite encounter time
    Fv_vc = logncdf(vc,mu,sigma);
    finiteTe(kk) = 1 - (1-Fv_vc^2)/4; 

    % Median encounter time -> Find zero of myfunc
    myfunc = @(t) computeCDFTe_Lognormal(mu,sigma,vc,dx0,t) - 0.5;

    % % Method 1 - Matlab's FZERO (unreliable for small sigV)
    % % Initial guess
    % % t0 = 2*dx0/muV;
    % % Median encounter time is zero of myfunc
    % % medianTe(kk) = fzero(myfunc,t0); 
    
    % Method 2 - BISECTION METHOD (slower but robust)
    % Upper and lower bound guess
    lb = 0.9*dx0/max(rangeVc);
    ub = 10*dx0/muV;
    
    a = lb; % by construction, myfunc(a) < 0
    b = ub; % by construction, myfunc(b) > 0

    keepgoing = true; 
    count = 0;
    while keepgoing
        c = (a+b)/2;
        if  myfunc(c) > 0; b = c; else; a = c; end

        count  = count + 1;
        relerr = abs(a-b)/(a+b);
        if count>1e5 || relerr < 1e-15
            % Stop search and save numerical solution
            keepgoing = false;
            medianTe(kk) = (a+b)/2;
        end
    end
end

if exist('outputfile','var')
    save(outputfile)
else
    save('output_autosave.mat','cdfTe', 'medianTe')
end
end