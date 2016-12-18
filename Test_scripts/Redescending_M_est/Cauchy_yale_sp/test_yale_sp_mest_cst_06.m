clear all;
close all;

% init

fname = 'cauchy_st';
noise_type = 'sp';
dataset = 'yale';

levels = [0.6];
opt.DEBUG = 1;

for n_level=levels

[O, X] = yale_sp(1, n_level);
X = X(:,:,1:64);
O = O(:,:,1:64);

[m, n, p] = size(X);
% lambda_opt = 1 / sqrt(max(m, n));
% binf = floor(log10(lambda_opt));
% bup = ceil(log10(lambda_opt));
% lambda = linspace(0.5*(10^binf), 1.5*(10^bup), 40) ;
lambdaS = min(size(X(:,:,1))) / sqrt(max(size(X(:,:,1))));
alpha = linspace(0, 1, 40);
sigma = logspace(log10(0.05), log10(2), 10);
H = ones(size(X));
X_cell = {zeros(size(X))};
eps = 1e-5;
max_iter = 500;

results_cauchy_st = cell(length(sigma), length(alpha));
i = 1;

for s=sigma
    j = 1;
    for a=alpha
%     opt.verbose = 1;
%     opt.lambda = l;
%     [L, S] = wrapper_georgios(@ten_rpca, X, opt);
        try
            tic
%             [D, I] = BRTF(X, k, nFrames, 'initVar', a, 'maxRank', 168);
            [L, iter,ranks] = lrtc_LiBCD_GS(X, H, X_cell, a*lambdaS,  s, 1, @grad_cauchy, @sv_shrinkage, eps, max_iter);
            time = toc;
            S = X - L;  % Heuristic

            msiqa = MSIQA(O, L);
            msiqa_f = MSIQA(O(:,:,1), L(:,:,1));

            results_cauchy_st{i,j} = struct(...
                'alpha', a, ...
                'sigma', s, ...
                'L', L, ...
                'E', S, ...
                'MSIQA', msiqa, ...
                'MSIQA_f', msiqa_f, ...
                'iter', 1:iter, ...
                'rank', ranks, ...
                'time', time ...
            );
        catch exception
            results_cauchy_st{i,j} = exception;
            fprintf('EXCEPTION: %s\n', exception.message);
        end

        j = j + 1;
    end
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_cauchy_st',...
         '-v7.3');

clear results_cauchy_st

end
