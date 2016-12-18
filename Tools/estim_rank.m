function [ r ] = estim_rank( A, t )
%ESTIM_RANK Estimates the rank of A based on the SVD
[~, S, ~] = svd(A);
S = diag(S);

if nargin < 2
    r = sum(S ./ max(S) > 1e-5);
else
    r = sum(S ./ max(S) > t);
end

end

