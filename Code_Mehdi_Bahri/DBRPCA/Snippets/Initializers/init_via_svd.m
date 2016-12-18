function [vars] = init_via_svd(vars, params)
%INIT_VIA_SVD ML estimation of the factors by SVD
%
% Mehdi Bahri - Imperial College London
% April, 2016

if params.TIME > 1
    tic
end

vars.T = zeros(params.r, params.r, params.Nobs);

% Left basis
vars.Uc = zeros(params.n, params.r);
% Right basis
vars.Ur = zeros(params.m, params.r);

for i=1:params.Nobs
    [U, S, V] = svd(vars.X(:,:,i));
    vars.Uc = vars.Uc + U(:,1:params.r);
    vars.Ur = vars.Ur + V(:,1:params.r);
    vars.T(:,:,i) = S(1:params.r, 1:params.r);
end

vars.Uc = vars.Uc / params.Nobs;
vars.Ur = vars.Ur / params.Nobs;

if params.TIME > 1
    svdtime = toc;
    fprintf('Initialized Uc, Ur, T via SVD (%fs)\n', svdtime);
end

end