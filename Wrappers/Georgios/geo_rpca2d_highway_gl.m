function [ Y, info ] = geo_rpca2d_highway_gl( algorithm, X, opt )
%GEO_RPCA2D_HIGHWAY Summary of this function goes here
%   Detailed explanation goes here

[Data, Info] = algorithm(X, 'lambda', ...
    opt.lambda, 'alpha_c', opt.p, 'alpha_r', opt.p, 'mean', true, 'tol', 1e-5);
Y.A = Data.L;
Y.E = Data.E;

info.err = Info.err;
info.iter = 1:Info.niter;
end

