function [ Data, Info ] = BRTF( X, k, nFrames, varargin )
%BRTF Calls BRTF and reconstructs the low-rank tensor

[imHeight, imWidth] = size(X(:,:,1));

Y = double(reshape(X, [imHeight*imWidth, k, nFrames]));

if max(Y(:)) > 1
    Y = rescale(Y, 0, 1);
end

[model] = BayesRCP(Y, 'maxiters', 500, varargin{:});

Data.E = model.E;
Data.L = double(ktensor(ones(model.TrueRank,1), model.Z{:}));

Data.L = squeeze(reshape(Data.L, [imHeight, imWidth, k, nFrames]));
Data.E = squeeze(reshape(Data.E, [imHeight, imWidth, k, nFrames]));

Info.LB = model.LowBound;
Info.TrueRank = model.TrueRank;
Info.gammas = model.gammas;
Info.err = model.err;

end

