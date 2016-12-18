clear all;
close all;

rng('default')

% init

fname = 'tnn';
noise_type = 'patch';
dataset = 'yale';

levels = [50, 100, 160];
opt.DEBUG = 1;

for n_level=levels

[O, X] = yale_patch(1, n_level);
X = X(:,:,1:64);
O = O(:,:,1:64);

[m, n, p] = size(X);
lambda_opt = 1 / sqrt(p * max(m, n));
binf = floor(log10(lambda_opt));
bup = ceil(log10(lambda_opt));
lambda = linspace(0.5*(10^binf), 1.5*(10^bup), 30) ;


results_trpca_tnn = cell(1, length(lambda));
i = 1;

for l=lambda
%     opt.verbose = 1;
%     opt.lambda = l;
%     [L, S] = wrapper_georgios(@ten_rpca, X, opt);
    tic
    [L, S, obj, err, iter] = trpca_tnn(X, l, opt);
    time = toc;

    msiqa = MSIQA(O, L);
    msiqa_f = MSIQA(O(:,:,1), L(:,:,1));

    results_trpca_tnn{i} = struct(...
        'lambda', l, ...
        'L', L, ...
        'E', S, ...
        'MSIQA', msiqa, ...
        'MSIQA_f', msiqa_f, ...
        'iter', iter, ...
        'err', err, ...
        'obj', obj, ...
        'time', time ...
    );
    
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_trpca_tnn',...
         '-v7.3');

clear results_trpca_tnn

end