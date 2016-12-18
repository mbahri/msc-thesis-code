function [X] = welsch_mest(X, sigma2)
% [X] = welsch_mest(X, sigma2)
% Computes the minimiser of the Welsch M-estimator for each element of X.
% INPUTS
%       X           an array of any dimensionality
%       sigma2      the square of the sigma parameter of the estimator
% OUTPUTS
%       X           the result
%
% Georgios Papamakarios
% Imperial College London
% July 2014

idx = X ~= 0;
X(idx) = X(idx) - X(idx) .* exp(- X(idx).^2 / sigma2);   
