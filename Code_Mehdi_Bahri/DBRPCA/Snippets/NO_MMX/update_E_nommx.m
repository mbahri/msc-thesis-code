function [ E ] = update_E_nommx(X, Y, T, Uc, Ur, params)
%UPDATE_E_NOMMX Solves for E by soft-thresholding
%
% Mehdi Bahri - Imperial College London
% April, 2016

E = zeros(params.n, params.m, params.Nobs);

mu = params.mu;
lambda = params.lambda;

if params.PARALLEL
    parfor k=1:params.Nobs
        E(:,:,k) = soft_shrinkage(X(:,:,k) - Uc*T(:,:,k)*Ur' + ...
            (1/mu)*Y(:,:,k), lambda / mu);
    end
else
    for k=1:params.Nobs
        E(:,:,k) = soft_shrinkage(X(:,:,k) - Uc*T(:,:,k)*Ur' + ...
            (1/mu)*Y(:,:,k), lambda / mu);
    end
end

end

