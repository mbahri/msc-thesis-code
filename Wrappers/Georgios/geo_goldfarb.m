function [ Y, info ] = geo_goldfarb( X, opt )
%GEO_GOLDFARB Summary of this function goes here
%   Detailed explanation goes here

[Data, Info] = rpca2d_l2(X, 'lambda', ...
    opt.lambda, 'mean', true);
Y.A = Data.L;
Y.E = Data.E;

info.err = Info.err;
info.iter = 1:Info.niter;
end

