clear all;
close all;

init

lambda = 2*logspace(-1, -4, 20);
alphas = logspace(0, -4, 5);
fname = 'mehdi_gl_l2';
noise_type = 'sp';
dataset = 'facade';

algorithms = {@rpca2d_gl_l2};
names = {'rpca2d_gl_l2'};

% Noise levels
levels = [0.1, 0.3, 0.6];

% Structure with all the information

for n_level = levels

%     [O, X] = yale_sp(1, n_level);
%     X = X(:,:,1:64);
%     O = O(:,:,1:64);

    [O, X] = facade_sp(1, n_level);
    if n_level < 0.3
        alpha_t = 1e-3;
    else
        alpha_t = 1e-2;
    end

    for a=1:length(algorithms);

        alg = algorithms{a};
        nam = names{a};
        
        results_mehdi.(nam) = cell(length(lambda), length(alphas));
        
        i = 1;
        for l=lambda
            j = 1;
            for alpha=alphas
                tic
                [Data, Info] = alg(X, ...
                    'lambda', l, ...
                    'mean', false, ...
                    'alpha_t', alpha_t, ...
                    'alpha_c', alpha, ...
                    'alpha_r', alpha, ...
                    'time', 0 ...       % Disable all timing informations
                );
                time = toc;

                L = Data.L;
                msiqa = MSIQA(O, L);
                msiqa_f = MSIQA(O(:,:,1), L(:,:,1));
                [fsim, fsimc] = FeatureSIM(O, L);

                results_mehdi.(nam){i, j} = struct(...
                    'lambda', l, ...
                    'alpha', alpha, ...
                    'L', L, ...
                    'E', Data.E, ...
                    'MSIQA', msiqa, ...
                    'fsim', fsim, ...
                    'fsimc', fsimc, ...
                    'MSIQA_f', msiqa_f, ...
                    'iter', 1:Info.niter, ...
                    'err', Info.err, ...
                    'time', time ...
                );
                j = j + 1;
            end
            i = i + 1;
        end
    end
    save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_mehdi',  '-v7.3');
    clear results_mehdi
end