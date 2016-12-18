function [ Y, info ] = geo_tensor_rpca( X, opt )
%GEO_TENSOR_RPCA Summary of this function goes here
%   Detailed explanation goes here

[L,S,err,iter] = tensor_rpca(X, opt.lambda);

Y.A = L;
Y.E = S;
info.iter = 1:iter;
info.err = err;


end

