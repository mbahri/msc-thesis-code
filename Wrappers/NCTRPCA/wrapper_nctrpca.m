function [ L, S ] = wrapper_nctrpca( X, r, t )
%WRAPPER_NCTRPCA Wrapper for non-convex tensor robust PCA
%
% Mehdi Bahri - Imperial College London
% July, 2016

% if nargin < 2
%     r = min(size(X(:,:,1)));
% end

params.thresh_const = t;

[L, S] = ncrpca_ten_new(X, r, params);

end