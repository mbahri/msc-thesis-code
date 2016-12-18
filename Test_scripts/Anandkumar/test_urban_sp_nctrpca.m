clear all;
close all;

% init

fname = 'nctrpca';
noise_type = 'sp';
dataset = 'urban';

levels = [0.1, 0.3, 0.6];
opt.DEBUG = 1;

for n_level=levels

[O, X] = urban_sp(0.5, n_level);

[m, n, p] = size(X);
% lambda_opt = 1 / sqrt(max(m, n));
% binf = floor(log10(lambda_opt));
% bup = ceil(log10(lambda_opt));
% lambda = linspace(0.5*(10^binf), 1.5*(10^bup), 40) ;

rankk = floor(linspace(25, 150, 10));
thr = logspace(1, 4, 4);


results = cell(length(thr), length(rankk));
i = 1;

for t=thr
    j = 1;
    for r=rankk
%     opt.verbose = 1;
%     opt.lambda = l;
%     [L, S] = wrapper_georgios(@ten_rpca, X, opt);
        try
            tic
%             [D, I] = BRTF(X, k, nFrames, 'initVar', a, 'maxRank', 168);
            [D, I] = wrapper_nctrpca(X, r, t);
            time = toc;

            msiqa = MSIQA(O, D.L);
            msiqa_f = MSIQA(O(:,:,1), D.L(:,:,1));

            results{i,j} = struct(...
                'rank', r, ...
                'threshold', t, ...
                'L', D.L, ...
                'E', D.E, ...
                'MSIQA', msiqa, ...
                'MSIQA_f', msiqa_f, ...
                'iter', 1:I.iter, ...
                'err', I.err, ...
                'time', time ...
            );
        catch exception
            results{i,j} = exception;
            fprintf('EXCEPTION: %s\n', exception.message);
        end

        j = j + 1;
    end
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results',...
         '-v7.3');

clear results

end