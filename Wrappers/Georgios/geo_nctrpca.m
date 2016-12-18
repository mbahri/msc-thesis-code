function [ Y, info ] = geo_nctrpca( X, opt )
%GEO_RPCA2D_HIGHWAY Summary of this function goes here
%   Detailed explanation goes here

[Data, Info] = wrapper_nctrpca(X, opt.lambda, opt.p);
Y.A = Data.L;
Y.E = Data.E;

info.err = Info.err;
info.iter = 1:Info.iter;
end

