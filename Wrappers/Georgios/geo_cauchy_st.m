function [ Y, info ] = geo_cauchy_st( X, opt )
%GEO_MEST Summary of this function goes here
%   Detailed explanation goes here

H = ones(size(X));
X_cell = {zeros(size(X))};

[L, iter,~] = lrtc_LiBCD_GS(X, H, X_cell, opt.lambda,  opt.p, ...
    1, @grad_cauchy, @sv_shrinkage, 1e-5, 500);

Y.A = L;
Y.E = X - L;
info.iter = 1:iter;
info.err = zeros(1, iter);

end

