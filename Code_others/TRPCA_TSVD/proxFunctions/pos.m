function P = pos(A)
% P = A .* double( A > 0 );
P = A;
P(P <= 0) = 0;
end