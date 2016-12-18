function [Y] = ten_shrinkage(X, t, p, pos)
% [Y] = ten_shrinkage(X, t, p, pos)
% Element-wise shrinkage operator on tensors.
% INPUTS
%       X      a tensor to be shrinked
%       t      the amount of shrinkage
%       p      the norm parameter for generalised shrinkage (optional)
%       pos    if true, the result can only be non-negative (optional)
% OUTPUTS
%       Y      the result of shrinking X
%
% Georgios Papamakarios
% Imperial College London
% May 2014

if nargin < 4
    pos = false;
end
if nargin < 3
    p = 1;
end

X = double(X);
Y = shrinkage(X, t, p, pos);
Y = tensor(Y);
