function [vars] = comp_Ur_mmx_L2(vars, params)
%COMP_UR_MMX_L2 Optimal Ur with squared Frob. penalty, MMX
%
% Mehdi Bahri - Imperial College London
% May, 2016

S = vars.S;
T = vars.T;
Uc = vars.Uc;

Urn = mmx('mult', S, Uc, 'tn');
Urn = sum(mmx('mult', Urn, T), 3);

V = params.mu*(Uc'*Uc);

Urd = mmx('mult', T, V, 'tn');
Urd = sum(mmx('mult', Urd, T), 3);

vars.Ur = Urn / (params.alpha_r * eye(params.r) + Urd);

end