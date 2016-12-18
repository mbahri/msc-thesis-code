function [ Y, Info ] = geo_rnndl( X, opt )
%GEO_RNNDL Wrapps NNDL for use with Georgio's benchmarks

params.alpha = opt.lambda;

if isfield(opt, 'thr')
    params.threshold = opt.thr
end

[h, w, N] = size(X);
XX = reshape(X, h*w, N);

[ D, Info ] = rnndl( XX, min([h, w, N]), params );
Y.A = reshape(D.R, h, w, N);
Y.E = X - Y.A;

end

