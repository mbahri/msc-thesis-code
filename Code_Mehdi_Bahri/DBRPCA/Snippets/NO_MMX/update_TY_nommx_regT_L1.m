function [T1, K1, Y1, Yt1] = update_TY_nommx_regT_L1(S, Uc, Ur, T, X, Y, Yt, E, params)
%UPDATE_TY_NOMMX_REGT_L1 Solves for T, L1 penalty, MMX
%
% Mehdi Bahri - Imperial College London
% July, 2016

tic

mu = params.mu;
mu_t = params.mu_t;
r = params.r;
alpha_t = params.alpha_t;

Y1 = zeros(size(X));
K1 = zeros(r, r, params.Nobs);

% K = (alpha_t * eye(r^2) + mu * kron(Ur'*Ur, Uc'*Uc));

U = (-mu / mu_t) * (Uc'*Uc);
V = Ur'*Ur;

if params.PARALLEL
    parfor k=1:params.Nobs
%         K1(:,:,k) = from_vec( K \ to_vec(Uc'*S(:,:,k)*Ur + mu_t * T + Yt) , r, r );
        K1(:,:,k) = dlyap(U, V, (Uc'*S(:,:,k)*Ur + mu_t * T(:,:,k) + Yt(:,:,k)) ./ mu_t );
        Y1(:,:,k) = Y(:,:,k) + mu*( X(:,:,k) - Uc*K1(:,:,k)*Ur' - E(:,:,k) );
    end
else
    for k=1:params.Nobs
%         K1(:,:,k) = from_vec( K \ to_vec(Uc'*S(:,:,k)*Ur + mu_t * T + Yt) , r, r );
        K1(:,:,k) = dlyap(U, V, (Uc'*S(:,:,k)*Ur + mu_t * T(:,:,k) + Yt(:,:,k)) ./ mu_t );
        Y1(:,:,k) = Y(:,:,k) + mu*( X(:,:,k) - Uc*K1(:,:,k)*Ur' - E(:,:,k) );
    end
end

vars.Y = Y1;
vars.K = K1;
vars.T = soft_shrinkage(K1 - Yt / mu_t, alpha_t/mu_t);
vars.Yt = Yt + mu_t*(T1 - K1);

up_t_time = toc;
if params.TIME
    fprintf('Updated L1 regularized T in %fs\n', up_t_time);
end

end

