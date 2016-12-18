function [arr] = mat2arr(mat, siz, rowidx, colidx)
% [arr] = mat2arr(mat, siz, rowidx, colidx)
% Transforms a matrix to an array.
% INPUTS
%       mat         the matrix to be transformed
%       siz         the size of the resulting array
%       rowidx      array indices that map to rows
%       colidx      array indices that map to columns (optional)
% OUTPUTS
%       arr         the output array
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

if nargin < 4
    colidx = setdiff(1:length(siz), rowidx);
end

mat = tenmat(mat, rowidx, colidx, siz);
arr = tensor(mat);
arr = double(arr);
