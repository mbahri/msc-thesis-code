function [ n ] = expectation_norm( Y, Uc, T, Ur, E, UcUc, UrUr, ...
    Sigma_T, S, Var_E )
%EXPECTATION_NORM Expectation of ||Yk - Uc*Tk*Ur' - Ek||_F^2

r = size(T, 1);

n = trace(Y'*Y) + 2*trace(E'*Y) - 2*trace(Y'*Uc*T*Ur') ...
    - 2*trace(E'*Uc*T*Ur');

n = n + trace(E'*E) + sum(Var_E(:));

% S = col2im(Sigma_T, [r, r], [r^2, r^2], 'distinct');
% S = blocktranspose(S, r^2, r^2, r, r);

% S = get_V(Sigma_T, r);

% n = n + to_vec(UcUc)' * ( kron(T, T) + S ) * to_vec(UrUr);
% n = n + vUcUc' * ( kron(T, T) + S ) * vUrUr;

n = n + quadsum(UcUc, UrUr, T, Sigma_T);

end

