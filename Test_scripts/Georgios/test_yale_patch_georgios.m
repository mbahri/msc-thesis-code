clear all;
close all;

init

lambda = 2*logspace(-1, -4, 30);
% lambda = 1e-1;
fname = 'georgios';
noise_type = 'patch';
dataset = 'yale';

algorithms = {@ten_rpca, @ten_rcpd, @ten_orpca, @ten_brpca};
names = {'ten_rpca', 'ten_rcpd', 'ten_orpca', 'ten_brpca'};

% Noise levels
levels = [50, 100, 160];

% Structure with all the information

for n_level = levels

    [O, X] = yale_patch(1, n_level);
    X = X(:,:,1:64);
    O = O(:,:,1:64);

%     facade_sp(1, n_level);

    for a=1:length(algorithms);

        alg = algorithms{a};
        nam = names{a};
        i = 1;
        
        results_georgios.(nam) = cell(1, length(lambda));

        for l=lambda
            % Parameters to set lambda and monitor the progress
            opt.verbose = 1;
            opt.lambda = l;
            
            % Retrieve the result and compute image fidelity measures
%             [L, S] = wrapper_georgios(alg, X, opt);
            tic
            [Y, info] = alg(X, opt);
            time = toc;
            
            L = double(Y.A);
            msiqa = MSIQA(O, L);
            msiqa_f = MSIQA(O(:,:,1), L(:,:,1));
            
            results_georgios.(nam){i} = struct(...
                'lambda', l, ...
                'L', L, ...
                'E', double(Y.E), ...
                'MSIQA', msiqa, ...
                'MSIQA_f', msiqa_f, ...
                'iter', info.iter, ...
                'err', info.err, ...
                'time', time ...
            );
            
            i = i + 1;
        end
    end
    
    save(sprintf('results_%s_%s_%s_%f.mat', dataset, fname, ...
        noise_type, n_level), ...
        'results_georgios',  '-v7.3');
    clear results_georgios
    
end