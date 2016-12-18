r = 10;

A = rand(r^2, r^2);
B = A'*A;

L = chol(B, 'lower');
% K = inv(L);
[~,K,~] = svd(B);

[X,Y] = krondecomp( K, r, r, r, r )

S = kron(X, Y);

norm(K - S, 'fro') / norm(K, 'fro')