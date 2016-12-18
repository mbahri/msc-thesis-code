function [n] = ew_norm(X, p) 
% [n] = ew_norm(X, p) 
% Element-wise Lp norm of array X.
% INPUTS
%       X     an array of any dimensionality
%       p     p parameter of the Lp norm
% OUTPUTS
%       n     element-wise Lp norm of X
%
% Georgios Papamakarios
% Imperial College London
% Apr 2014

X = X(:);
n = norm(X, p);
