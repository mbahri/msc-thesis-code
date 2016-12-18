function [T1, Y1] = update_TY_nommx_regT_L2(S, Uc, Ur, X, Y, E, params)
%UPDATE_TY_NOMMX_REGT_L2 Solves for T with a squared Frobenius penalty
%
% Mehdi Bahri - Imperial College London
% July, 2016

tic

mu = params.mu;
r = params.r;
alpha_t = params.alpha_t;
T1 = zeros(r, r, params.Nobs);
Y1 = zeros(size(X));

% K = (alpha_t * eye(r^2) + mu * kron(Ur'*Ur, Uc'*Uc));

U = (-mu/alpha_t)*(Uc'*Uc);
V = Ur'*Ur;

if params.PARALLEL
    parfor k=1:params.Nobs
%         T1(:,:,k) = from_vec( K \ to_vec(Uc'*S(:,:,k)*Ur) , r, r );
        T1(:,:,k) = dlyap(U, V, (Uc'*S(:,:,k)*Ur) ./ alpha_t);
        Y1(:,:,k) = Y(:,:,k) + mu*( X(:,:,k) - Uc*T1(:,:,k)*Ur' - E(:,:,k) );
    end
else
    for k=1:params.Nobs
%         T1(:,:,k) = from_vec( K \ to_vec(Uc'*S(:,:,k)*Ur) , r, r );
        T1(:,:,k) = dlyap(U, V, (Uc'*S(:,:,k)*Ur) ./ alpha_t);
        Y1(:,:,k) = Y(:,:,k) + mu*( X(:,:,k) - Uc*T1(:,:,k)*Ur' - E(:,:,k) );
    end
end

up_t_time = toc;
if params.TIME
    fprintf('Updated L2 regularized T in %fs\n', up_t_time);
end

end

