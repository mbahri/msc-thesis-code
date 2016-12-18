function [X] = l1_l2_mest(X, alpha, sigma2)
% [X] = l1_l2_mest(X, alpha, sigma2)
% Computes the minimiser of the L1-L2 M-estimator for each element of X.
% INPUTS
%       X           an array of any dimensionality
%       alpha       the alpha parameter of the estimator
%       sigma2      the square of the sigma parameter of the estimator
% OUTPUTS
%       X           the result
%
% Georgios Papamakarios
% Imperial College London
% July 2014

idx = X ~= 0;
X(idx) = X(idx) - X(idx) ./ sqrt(alpha + X(idx).^2 / sigma2);
