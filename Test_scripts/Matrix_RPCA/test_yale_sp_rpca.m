clear all;
close all;

% init

fname = 'rpca';
noise_type = 'sp';
dataset = 'yale';

levels = [0.1, 0.3, 0.6];
opt.DEBUG = 1;

for n_level=levels

[O, X] = yale_sp(1, n_level);
X = X(:,:,1:64);
O = O(:,:,1:64);

[m, n, p] = size(X);
lambda_opt = 1 / sqrt(max(m, n));
binf = floor(log10(lambda_opt))-2;
bup = ceil(log10(lambda_opt));
lambda = linspace(0.5*(10^binf), 1.5*(10^bup), 100) ;


results_rpca = cell(1, length(lambda));
i = 1;

for l=lambda
    tic
    [D, I] = wrapper_rpca(X, l);
    time = toc;

    msiqa = MSIQA(O, D.L);
    msiqa_f = MSIQA(O(:,:,1), D.L(:,:,1));
    [fsim, fsimc] = FeatureSIM(O, D.L);

    results_rpca{i} = struct(...
        'lambda', l, ...
        'L', D.L, ...
        'E', D.E, ...
        'MSIQA', msiqa, ...
        'fsim', fsim, ...
        'fsimc', fsimc, ...
        'MSIQA_f', msiqa_f, ...
        'iter', 1:I.iter, ...
        'rank', I.rank ...
    );
    
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_rpca',...
         '-v7.3');

clear results_rpca

end