function [ err ] = matrix_relative_error( A, B )
%MATRIX_RELATIVE_ERROR Relative error between A and B in Frob. norm

err = norm(A - B, 'fro') / norm(A, 'fro');

end

