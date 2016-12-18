clear;
close all;
% clc;

synthdata_2;

lambda = [1 8:1e-1:12 20 50];
mu = 5:5:110;

% lambda = [10];
% mu = 1e3;

Nl = length(lambda);
Nm = length(mu);

E_rec = zeros(length(mu), length(lambda));
O_rec = zeros(length(mu), length(lambda));
Full_rec = zeros(length(mu), length(lambda));

stemUc = cell(length(mu), length(lambda));
stemUr = cell(length(mu), length(lambda));
% meshUc = cell(length(mu), length(lambda));
% meshUr = cell(length(mu), length(lambda));
% stemT = cell(length(mu), length(lambda));

for ll=1:Nl
    for mm=1:Nm
        l = lambda(ll);
        m = mu(mm);
        
        [ Uc, Ur, T, E ] = brpca( X, 20, 'rho', 2.9, 'maxiter', 1000, ...
            'mu1', m, 'mu2', m, 'lambda', l);
        stemUc{mm, ll} = stemplot(Uc, m, l);
        stemUr{mm, ll} = stemplot(Ur, m, l);
%         meshUc{mm, ll} = meshplot(Uc, m, l);
%         meshUr{mm, ll} = meshplot(Ur, m, l);
%         stemT{mm, ll} = stemplot(T(:,:,1), m, l);
        
        E_rec(mm, ll) = max_rec_e(E, E_ori);
        O_rec(mm, ll) = max_rec_o(Uc, T, Ur, O);
        Full_rec(mm, ll) = max_rec_f(Uc, T, Ur, E, X);
    end
end
        
save('dump_all_data_fine')