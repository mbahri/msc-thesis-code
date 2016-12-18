clear all;
close all;

% init

fname = 'tsvd';
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
binf = floor(log10(lambda_opt));
bup = ceil(log10(lambda_opt));
lambda = linspace(0.5*(10^binf), 1.5*(10^bup), 40) ;


results_trpca_tsvd = cell(1, length(lambda));
i = 1;

for l=lambda
%     opt.verbose = 1;
%     opt.lambda = l;
%     [L, S] = wrapper_georgios(@ten_rpca, X, opt);
    try
        tic
        [L, S, err, iter] = tensor_rpca(X, l);
        time = toc;

        msiqa = MSIQA(O, L);
        msiqa_f = MSIQA(O(:,:,1), L(:,:,1));

        results_trpca_tsvd{i} = struct(...
            'lambda', l, ...
            'L', L, ...
            'E', S, ...
            'MSIQA', msiqa, ...
            'MSIQA_f', msiqa_f, ...
            'iter', iter, ...
            'err', err, ...
            'time', time ...
        );
    catch exception
        results_trpca_tsvd{i} = exception;
        fprintf('EXCEPTION: %s\n', exception.message);
    end
    
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_trpca_tsvd',...
         '-v7.3');

clear results_trpca_tsvd

end