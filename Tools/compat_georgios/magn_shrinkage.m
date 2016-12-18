function [Y] = magn_shrinkage(X, t)
% [Y] = magn_shrinkage(X, t)
% Magnitude shrinkage operator on arrays, based on the elementwise L2 norm.
% INPUTS
%       X      an array of any dimensionality
%       t      the amount of shrinkage
% OUTPUTS
%       Y      the result of shrinking X
%
% Georgios Papamakarios
% Imperial College London
% June 2014

normX = ew_norm(X, 2);

if normX > t
    Y = X * (normX - t) / normX;
else
    Y = zeros(size(X));
end
