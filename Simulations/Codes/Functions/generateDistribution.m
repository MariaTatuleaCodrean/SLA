function [vel_pdf,mu,sigma] = generateDistribution(muV,sigV)
% GENERATEDISTRIBUTION generates a Lognormal distribution using the desired 
% mean and variance of the speed distribution.

% Calculate Lognormal distribution parameters from muV and sigmaV
sig2V   = sigV^2;
mu      = log(muV^2./sqrt(muV^2 + sig2V));
sig2    = log(1+sig2V/muV^2);
sigma   = sqrt(sig2);

% Generate lognormal distribution with those parameters
vel_pdf = makedist('Lognormal','mu',mu,'sigma',sigma);
end