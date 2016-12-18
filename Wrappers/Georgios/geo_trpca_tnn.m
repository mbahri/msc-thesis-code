function [ Y, info ] = geo_trpca_tnn( X, opt )
%GEO_TRPCA_TNN Summary of this function goes here
%   Detailed explanation goes here

opts.DEBUG = 1;

[L,S,obj,err,iter] = trpca_tnn(X, opt.lambda, opts);

Y.A = L;
Y.E = S;
info.iter = 1:iter;
info.err = err;

end

