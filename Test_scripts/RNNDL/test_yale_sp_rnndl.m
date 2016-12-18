clear all;
close all;

% init

fname = 'rnndl';
noise_type = 'sp';
dataset = 'yale';

levels = [0.1, 0.3, 0.6];
opt.DEBUG = 1;

for n_level=levels

[O, X] = yale_sp(1, n_level);
X = X(:,:,1:64);
O = O(:,:,1:64);

[m, n, p] = size(X);

lambda = 0:0.5:5;


results_rnndl = cell(1, length(lambda));
i = 1;

for l=lambda
    tic
    opt.lambda = l;
    [D, I] = geo_rnndl(X, opt);
    time = toc;

    msiqa = MSIQA(O, D.A);
    msiqa_f = MSIQA(O(:,:,1), D.A(:,:,1));
    [fsim, fsimc] = FeatureSIM(O, D.A);

    results_rnndl{i} = struct(...
        'lambda', l, ...
        'L', D.A, ...
        'E', D.E, ...
        'MSIQA', msiqa, ...
        'fsim', fsim, ...
        'fsimc', fsimc, ...
        'MSIQA_f', msiqa_f, ...
        'iter', I.iter ...
    );
    
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_rnndl',...
         '-v7.3');

clear results_rnndl

end