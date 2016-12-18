function [ Y, info ] = geo_rpca2d_yale( X, opt )
%GEO_RPCA2D_HIGHWAY Summary of this function goes here
%   Detailed explanation goes here

[Data, Info] = rpca2d_l2(X, 'lambda', ...
    opt.lambda, 'mean', false);
Y.A = Data.L;
Y.E = Data.E;

info.err = Info.err;
info.iter = 1:Info.niter;
end

