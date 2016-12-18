function [ vars ] = comp_Ur_nommx_L2( vars, params )
%COMP_UR_NOMMX_L2 Computes the optimal Ur with squared Frob. penalty
%
% Mehdi Bahri - Imperial College London
% April, 2016

S = vars.S;
T = vars.T;
Uc = vars.Uc;

mu = params.mu;

Urn = zeros(params.m, params.r);
Urd = zeros(params.r, params.r);

if params.PARALLEL
    parfor k=1:params.Nobs
        Urn = Urn + S(:,:,k)' * Uc * T(:,:,k);
        Urd = Urd + mu*T(:,:,k)'*(Uc'*Uc)*T(:,:,k);
    end
else
    for k=1:params.Nobs
        Urn = Urn + S(:,:,k)' * Uc * T(:,:,k);
        Urd = Urd + mu*T(:,:,k)'*(Uc'*Uc)*T(:,:,k);
    end
end

vars.Ur = Urn / (params.alpha_r*eye(params.r) + Urd);

end

