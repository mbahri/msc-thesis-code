function [ Y ] = blocktranspose( X, m, n, p, q )
%BLOCKTRANSPOSE Block transpose a matrix X

Y = reshape(X, [ p m/p q n/q ]);
Y = permute(Y, [ 1 4 3 2 ]);
Y = reshape(Y, [ p*n/q q*m/p]);

end

