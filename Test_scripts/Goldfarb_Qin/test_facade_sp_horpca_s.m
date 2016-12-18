% clear all;
% close all;

% init

fname = 'horpca_s';
noise_type = 'sp';
dataset = 'facade';

levels = [0.1, 0.3, 0.6];
opt.DEBUG = 1;

for n_level=levels

[O, X] = facade_sp(1, n_level);

[m, n, p] = size(X);
lambda_opt = 1 / sqrt(max(size(X)));
lambda = logspace(-1, 1, 30) * lambda_opt;


results_horpca_s.horpca_s = cell(1, length(lambda));
results_horpca_s.horpca_s_tc = cell(1, length(lambda));
i = 1;

for l=lambda
%     opt.verbose = 1;
%     opt.lambda = l;
%     [L, S] = wrapper_georgios(@ten_rpca, X, opt);

    % First the normal algorithm
    try
        tic
        [D, I] = horpca_s(X, l);
        time = toc;

        msiqa = MSIQA(O, D.L);
        msiqa_f = MSIQA(O(:,:,1), D.L(:,:,1));

        results_horpca_s.horpca_s{i} = struct(...
            'lambda', l, ...
            'L', D.L, ...
            'E', D.E, ...
            'MSIQA', msiqa, ...
            'MSIQA_f', msiqa_f, ...
            'iter', 1:I.iter, ...
            'err', I.err, ...
            'time', time ...
        );
    catch exception
        results_horpca_s.horpca_s{i} = exception;
        fprintf('EXCEPTION: %s\n', exception.message);
    end
    
    % Then with missing values enabled
     try
        tic
        [D, I] = horpca_s(X, l, true);
        time = toc;

        msiqa = MSIQA(O, D.L);
        msiqa_f = MSIQA(O(:,:,1), D.L(:,:,1));

        results_horpca_s.horpca_s_tc{i} = struct(...
            'lambda', l, ...
            'L', D.L, ...
            'E', D.E, ...
            'MSIQA', msiqa, ...
            'MSIQA_f', msiqa_f, ...
            'iter', 1:I.iter, ...
            'err', I.err, ...
            'time', time ...
        );
    catch exception
        results_horpca_s.horpca_s_tc{i} = exception;
        fprintf('EXCEPTION: %s\n', exception.message);
    end
    
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), 'results_horpca_s',  '-v7.3');

clear results_horpca_s

end