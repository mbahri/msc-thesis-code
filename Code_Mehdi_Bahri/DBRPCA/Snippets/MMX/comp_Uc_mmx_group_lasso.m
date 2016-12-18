function [vars] = comp_Uc_mmx_group_lasso(vars, params)
%COMP_UC_MMX_GROUP_LASSO Optimal Uc with Group Lasso penalty, MMX
%
% Mehdi Bahri - Imperial College London
% June, 2016

S = vars.S;
T = vars.T;
Ur = vars.Ur;
Y_c = vars.Y_c;

Ucn = mmx('mult', S, Ur);
Ucn = sum(mmx('mult', Ucn, T, 'nt') ,3) + vars.mu_c*vars.A + Y_c;

V = params.mu*(Ur'*Ur);

Ucd = mmx('mult', T, V);
Ucd = sum(mmx('mult', Ucd, T, 'nt'), 3);

Uc = Ucn / (vars.mu_c * eye(params.r) + Ucd);
vars.Uc = Uc;
vars.A = solve_l1l2(Uc - Y_c ./ vars.mu_c, params.alpha_c / vars.mu_c);

vars.Y_c = Y_c + vars.mu_c * (vars.A - Uc); % Need to check the sign
vars.mu_c = params.rho*vars.mu_c;

end