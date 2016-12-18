% clear

init

% synthdata_2

% load_X
% Copy_2_of_load_X
% E_ori = zeros(128,128,2414);

r = 21;
maxiter = 1;

% [Ucd, Urd, Td, Ed] = brpca_mmx(X, r, 'alpha_c', 1, 'alpha_r', 1, 'maxiter', 1000, 'rho', 11);

% Ucd = rand(n,r);
% Urd = rand(m,r);
% Td = rand(r,r,Nobs);
% Ed = rand(n,m,Nobs);

estim_rank(Ucd)
estim_rank(Urd)
max_rec_e(Ed, E_ori)
max_rec_o(Ucd, Td, Urd, O)
max_rec_f(Ucd, Td, Urd, Ed, X)

% stemplot(Ucd,1,1);
% waitforbuttonpress;
% stemplot(Urd,1,1);
% waitforbuttonpress;
% meshplot(Ucd,1,1);
% waitforbuttonpress;
% meshplot(Urd,1,1);
% waitforbuttonpress;
% meshplot(Td(:,:,1),1,1);
% waitforbuttonpress;
% meshplot(Ed(:,:,1),1,1);
% waitforbuttonpress;
% fprintf('Continuing\n');

% [ Ucv1, Urv1, Tv1, Ev1, PUcv1, PUrv1 ] = vbbrpca_kronfac_nommx( X, Td, Ed, Ucd, Urd, r, maxiter, O, E_ori);
% [ Ucv2, Urv2, Tv2, Ev2, PUcv2, PUrv2 ] = vbbrpca_mmx( X, Td, Ed, Ucd, Urd, r, maxiter, O, E_ori);
[ Ucv3, Urv3, Tv3, Ev3, PUcv3, PUrv3 ] = vbbrpca( X, Td, Ed, Ucd, Urd, r, maxiter, O, E_ori);
% [ Ucv3, Urv3, Tv3, Ev3, PUcv3, PUrv3 ] = vbbrpca_p( X, Td, Ed, Ucd, Urd, r, maxiter, O, E_ori);

% [ Ucv2, Urv2, Tv2, Ev2, PUcv2, PUrv2 ] = vbbrpca_mmx( X(:,:,1), Td(:,:,1),...
%     Ed(:,:,1), Ucd, Urd, r, maxiter, O(:,:,1), E_ori(:,:,1));
