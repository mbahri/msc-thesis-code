function [ Y, info ] = geo_brtf_hall( X, opt )
%GEO_BRTF Summary of this function goes here
%   Detailed explanation goes here

[ Data, Info ] = BRTF( X, 1, size(X,3), 'initVar', opt.lambda, 'updateHyper', 'off', 'maxRank', 100 );
Y.A = Data.L;
Y.E = Data.E;

info.err = Info.err;
info.iter = 1:length(Info.LB);

end

