function [vars] = update_TY_mmx_regT_L2(vars, params) %S, Uc, Ur, X, Y, E
%UPDATE_TY_MMX_REGT_L2 Solves for T. Squared Frobenius penalty, MMX
%
% Mehdi Bahri - Imperial College London
% July, 2016

if params.TIME > 2
    tic
end

mu = params.mu;
r = params.r;
alpha_t = params.alpha_t;
Uc = vars.Uc;
Ur = vars.Ur;
Xt = vars.Xt;    % Use X tilde
S = vars.S;
Y = vars.Y;

T1 = zeros(r, r, params.Nobs);

% K = (alpha_t * eye(r^2) + mu * kron(Ur'*Ur, Uc'*Uc));
red_S = mmx('mult', Uc, S, 'tn');
% red_S = mmx('mult', red_S, Ur);
red_S = mmx('mult', red_S, Ur) ./ alpha_t;

U = (-mu/alpha_t)*(Uc'*Uc);
V = Ur'*Ur;

for k=1:params.Nobs
    T1(:,:,k) = dlyap(U, V, red_S(:,:,k));
end

% L = double(tenmat(tensor(red_S), 3))';
% temp_T = K \ L;
% for k=1:params.Nobs
%     T1(:,:,k) = from_vec(temp_T(:,k), r, r);
% end

A = mmx('mult', Uc, T1);
A = mmx('mult', A, Ur, 'nt');

% Update the structure
vars.Y = Y + mu*(Xt - A);
vars.T = T1;

if params.TIME > 2
    up_t_time = toc;
    fprintf('Updated L2 regularized T in %fs (MMX)\n', up_t_time);
end

end

