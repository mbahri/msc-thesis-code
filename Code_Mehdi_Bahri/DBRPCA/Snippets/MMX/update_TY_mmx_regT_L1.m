function [vars] = update_TY_mmx_regT_L1(vars, params)
%UPDATE_TY_MMX_REGT_L1 Solves for T, L1 penalty, MMX
%
% Mehdi Bahri - Imperial College London
% July, 2016

if params.TIME > 2
    tic
end

S = vars.S;
Uc = vars.Uc;
Ur = vars.Ur;
Xt = vars.Xt;
Y = vars.Y;
Yt = vars.Yt;

mu = params.mu;
mu_t = vars.mu_t;
r = params.r;
alpha_t = params.alpha_t;

T1 = zeros(r, r, params.Nobs);

red_S = mmx('mult', Uc, S, 'tn');
red_S = ( mmx('mult', red_S, Ur) + mu_t * vars.K + Yt ) ./ mu_t;

U = (-mu / mu_t) * (Uc'*Uc);
V = Ur'*Ur;

for k=1:params.Nobs
   T1(:,:,k) = dlyap(U, V, red_S(:,:,k));
end

K1 = soft_shrinkage(T1 - Yt / mu_t, alpha_t/mu_t);

vars.T = T1;
vars.K = K1;

A = mmx('mult', Uc, T1);
A = mmx('mult', A, Ur, 'nt');

vars.Y = Y + mu*(Xt - A);
vars.Yt = Yt + mu_t*(K1 - T1);
vars.mu_t = params.rho * vars.mu_t;

if params.TIME > 2
    up_t_time = toc;
    fprintf('Updated L1 regularized T in %fs (MMX)\n', up_t_time);
end

end

