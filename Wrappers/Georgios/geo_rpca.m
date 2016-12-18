function [ Data, Info ] = geo_rpca( X, opt )
%GEO_RPCA Wrapps the code for RPCA with Inexact ALM
warning('off', 'PROPACK:NotUsingMex');

[m, n, k] = size(X);
D = reshape(X, [m*n, k, 1]);

if nargin < 2
    [~, A_hat, E_hat, iter] = rpca(D);    
else
    [~, A_hat, E_hat, iter] = rpca(D, opt.lambda);
end

Info.iter = 1:iter;
Info.err = zeros(1, iter);

Data.A = reshape(A_hat, [m, n, k, 1]);
Data.E = reshape(E_hat, [m, n, k, 1]);

end

