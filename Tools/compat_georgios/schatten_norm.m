function [n] = schatten_norm(X, p)
% [n] = schatten_norm(X, p)
% Schatten p-norm of matrix X.
% INPUTS
%       X     a matrix
%       p     p parameter of the Schatten p-norm
% OUTPUTS
%       n     Schatten p-norm of X
%
% Georgios Papamakarios
% Imperial College London
% June 2014

[n1, n2] = size(X);

if n1 > n2
    Y = X' * X;
else
    Y = X * X';
end

n = trace(Y ^ (p/2));
n = n ^ (1/p);
