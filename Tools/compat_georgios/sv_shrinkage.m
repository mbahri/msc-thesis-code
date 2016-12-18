function [Y] = sv_shrinkage(X, t, p)
% [Y] = sv_shrinkage(X, t, p)
% Singular value shrinkage operator on matrices.
% INPUTS
%       X      a real matrix
%       t      the amount of shrinkage
%       p      the norm parameter for generalised shrinkage (optional)
% OUTPUTS
%       Y      the result of shrinking X
%
% Georgios Papamakarios
% Imperial College London
% Apr 2014

if all(t == 0)
    
    Y = X;

else
    
    if nargin < 3
        p = 1;
    end

    X = double(X);
    [U, S, V] = svd(X, 'econ');
    S = diag(shrinkage(diag(S), t, p, true));
    
    Y = U * S * V';
end