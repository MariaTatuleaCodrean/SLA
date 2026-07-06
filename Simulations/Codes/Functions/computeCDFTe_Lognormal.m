function [output] = computeCDFTe_Lognormal(mu,sigma,vc,dx0,t)
%COMPUTECDFTE_LOGNORMAL Cumulative distribution function of encounter time
%   Given a lognormal distribution (mu,sigma) and velocity of the central
%   bacterium (vc), this function calculates the probability that the first
%   encounter time to either its front or rear neighbour, at given distance
%   apart (dx0), is smaller than some time t (elements of t vector)

    % Intermediate steps
    Fv_plus  = logncdf((vc + dx0./t),mu,sigma);
    Fv_minus = logncdf(abs(vc - dx0./t),mu,sigma).*sign(vc - dx0./t);

    % P(rear encounter time < t)
    F_Tr = (1-Fv_plus)/2;
    % P(front encounter time < t)
    F_Tf = (1 + Fv_minus)/2;
        
    % Evaluate CDF of encounter time 
    output =  F_Tr + F_Tf - F_Tr.*F_Tf;

end