function [ v ] = to_vec( M )
%TO_VEC Implements the vec operator
[n, m] = size(M);

v = reshape(M, n*m, 1);

end

