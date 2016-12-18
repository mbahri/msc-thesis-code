function [vars] = comp_Uc_mmx_L2(vars, params) %S, T, Ur, 
%COMP_UC_MMX_L2 Optimal Uc with squared Frob. penalty
%
% Mehdi Bahri - Imperial College London
% May, 2016

S = vars.S;
T = vars.T;
Ur = vars.Ur;

Ucn = mmx('mult', S, Ur);
Ucn = sum(mmx('mult', Ucn, T, 'nt') ,3);

V = params.mu*(Ur'*Ur);

Ucd = mmx('mult', T, V);
Ucd = sum(mmx('mult', Ucd, T, 'nt'), 3);

vars.Uc = Ucn / (params.alpha_c * eye(params.r) + Ucd);

end