function [mat] = arr2mat(arr, rowidx, colidx)
% [mat] = arr2mat(arr, rowidx, colidx)
% Transforms an array to a matrix.
% INPUTS
%       arr         the array to be transformed
%       rowidx      array indices that map to rows
%       colidx      array indices that map to columns (optional)
% OUTPUTS
%       mat         the output matrix
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

if nargin < 3
    colidx = setdiff(1:length(size(arr)), rowidx);
end

arr = tensor(arr);
mat = tenmat(arr, rowidx, colidx);
mat = double(mat);
