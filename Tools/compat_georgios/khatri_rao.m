function [Y] = khatri_rao(X)
% [Y] = khatri_rao(X)
% Calculates the Khatri-Rao product of a set of matrices.
% INPUTS
%       X     cell array with matrices of the same number of columns
% OUTPUTS
%       Y     the Khatri-Rao product of the matrices in X
%
% Georgios Papamakarios
% Imperial College London
% May 2014

Y = X{1};

for i = 2:length(X)
    Y = khatri_rao_once(Y, X{i});
end


function [Y] = khatri_rao_once(A, B)
% Calculates the Khatri-Rao product of two matrices.
% INPUTS
%       A,B   matrices with the same number of columns
% OUTPUTS
%       Y     the Khatri-Rao product of A and B

N = size(A,2);

if N ~= size(B,2)
    error('Matrices must have the same number of columns.');
end

Y = zeros(size(A,1) * size(B,1), N);

for i = 1:N
    temp = B(:,i) * A(:,i)';
    Y(:,i) = temp(:);
end
