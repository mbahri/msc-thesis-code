function [ Y, info ] = geo_horpca( X, opt )
%GEO_HORPCA Summary of this function goes here
%   Detailed explanation goes here

[Data, Info] = horpca_s(X, opt.lambda);

Y.A = Data.L;
Y.E = Data.E;

info.err = Info.err;
info.iter = Info.iter;

end

