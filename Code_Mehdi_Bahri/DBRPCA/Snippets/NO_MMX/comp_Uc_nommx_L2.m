function [ vars ] = comp_Uc_nommx_L2( vars, params )
%COMP_UC_NOMMX_L2 Computes the optimal Ur with squared Frob. penalty
%
% Mehdi Bahri - Imperial College London
% April, 2016

S = vars.S;
T = vars.T;
Ur = vars.Ur;

mu = params.mu;

Ucn = zeros(params.n, params.r);
Ucd = zeros(params.r, params.r);

if params.PARALLEL
    parfor k=1:params.Nobs
        Ucn = Ucn + S(:,:,k) * Ur * T(:,:,k)';
        Ucd = Ucd + mu*T(:,:,k)*(Ur'*Ur)*T(:,:,k)';
    end
else
    for k=1:params.Nobs
        Ucn = Ucn + S(:,:,k) * Ur * T(:,:,k)';
        Ucd = Ucd + mu*T(:,:,k)*(Ur'*Ur)*T(:,:,k)';
    end
end
    
vars.Uc = Ucn / (params.alpha_c*eye(params.r) + Ucd);

end

