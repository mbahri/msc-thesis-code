clear all;
close all;

% init

fname = 'brtf';
noise_type = 'sp';
dataset = 'facade';
k = 3;
nFrames = 1;

levels = [0.1, 0.3, 0.6];
opt.DEBUG = 1;

for n_level=levels

[O, X] = facade_sp(1, n_level);

[m, n, p] = size(X);
% lambda_opt = 1 / sqrt(max(m, n));
% binf = floor(log10(lambda_opt));
% bup = ceil(log10(lambda_opt));
% lambda = linspace(0.5*(10^binf), 1.5*(10^bup), 40) ;
lambda = logspace(-3, 3, 20);   % actually initVar


results_brtf = cell(1, length(lambda));
i = 1;

for l=lambda
%     opt.verbose = 1;
%     opt.lambda = l;
%     [L, S] = wrapper_georgios(@ten_rpca, X, opt);
    try
        tic
        [D, I] = BRTF(X, k, nFrames, 'initVar', l, 'maxRank', 168);
        time = toc;
        L = D.L;
        S = D.E;

        msiqa = MSIQA(O, L);
        msiqa_f = MSIQA(O(:,:,1), L(:,:,1));

        results_brtf{i} = struct(...
            'lambda', l, ...
            'L', L, ...
            'E', S, ...
            'MSIQA', msiqa, ...
            'MSIQA_f', msiqa_f, ...
            'iter', 1:length(I.LB), ...
            'err', I.err, ...
            'gammas', I.gammas, ...
            'rank', I.TrueRank, ...
            'time', time ...
        );
    catch exception
        results_brtf{i} = exception;
        fprintf('EXCEPTION: %s\n', exception.message);
    end
    
    i = i + 1;
end

save( sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_brtf',...
         '-v7.3');

clear results_brtf

end