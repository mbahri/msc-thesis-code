function [ L, S ] = wrapper_tensor_rpca( X, lambda )
%WRAPPER_TENSOR_RPCA Calls the code for the method 'Novel methods for image
%completion and denoising based on tensor SVD'
%
% Mehdi Bahri - Imperial College London
% July, 2016

if nargin < 2
    [m, n, ~] = size(X);
    lambda = 1 / sqrt(max(m, n));
    fprintf('Lambda is %g\n', lambda);
    [L, S] = tensor_rpca(X, lambda);
else
    fprintf('Lambda is %g\n', lambda);
    [L, S] = tensor_rpca(X, lambda);
end

end

