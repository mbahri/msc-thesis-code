function [ L, S ] = wrapper_trpca_tnn( X, lambda, opts )
%WRAPPER_TRPCA_TNN Calls...
%
% Mehdi Bahri - Imperial College London
% July, 2016

if nargin < 2
    [m, n, p] = size(X);
    lambda = 1 / sqrt(p * max(m, n));
    OO.DEBUG = 1;
    fprintf('Lambda is %g\n', lambda);
    [L, S] = trpca_tnn(X, lambda, OO);
elseif nargin < 3
    OO.DEBUG = 1;
    fprintf('Lambda is %g\n', lambda);
    [L, S] = trpca_tnn(X, lambda, OO);
else
   OO = opts;
   if ~isfield(OO, 'DEBUG')
       OO.DEBUG = 1;
   end
   fprintf('Lambda is %g\n', lambda);
   [L, S] = trpca_tnn(X, lambda, OO);
end

end

