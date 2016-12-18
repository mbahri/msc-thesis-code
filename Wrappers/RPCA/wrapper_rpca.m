function [ Data, Info ] = wrapper_rpca( X, lambda )
%WRAPPER_RPCA Wrapps the code for RPCA with Inexact ALM
warning('off', 'PROPACK:NotUsingMex');

[m, n, k] = size(X);
D = reshape(X, [m*n, k, 1]);

if nargin < 2
    [EstimatedRank, A_hat, E_hat, iter] = rpca(D);    
else
    [EstimatedRank, A_hat, E_hat, iter] = rpca(D, lambda);
end

Info.rank = EstimatedRank;
Info.iter = iter;

Data.L = reshape(A_hat, [m, n, k, 1]);
Data.E = reshape(E_hat, [m, n, k, 1]);
Data.A_h = A_hat;
Data.E_h = E_hat;

end

