% experiment_highway_no_noise.m
% Executes the background subtraction experiment on the highway video,
% without the addition of noise.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% directory where to write
rootdir = '/vol/bitbucket/mb2215/experiments/highway';
matdir = fullfile(rootdir, 'matrices');
tendir = fullfile(rootdir, 'tensors');
mkdir(matdir);
mkdir(tendir);

% data parameters
datapath = fullfile('..', 'data');
datafile = 'highway.mat';
numims = 400;
% load(fullfile(datapath, datafile));

% MODIFICATION BY MEHDI BAHRI
load_highway
rowidx = [1, 2];
colidx = 3;

% algorithm parameters
lambda = logspace(-1, -4, 7);
p = 1;
mat_rank = 10 : 10 : 50;
ten_rank = 10 : 10 : 50;
ten_nrank_perc = 0.2 : 0.2 : 1.0;
ten_nrank = cell(1, length(ten_nrank_perc));

for ni = 1:length(numims)
    
    thisX = X(:, :, 1:numims(ni));
        
    for i = 1:length(ten_nrank_perc)
        ten_nrank{i} = floor(ten_nrank_perc(i) * size(thisX));
    end
    
    mat_data.A = arr2mat(thisX, rowidx, colidx);
    ten_data.A = thisX;
    
    % no noise
    mat_data.X = mat_data.A;
    ten_data.X = ten_data.A;

    filename = sprintf('%dx%dx%d_no_noise_part2_brtf', size(thisX));

%     fprintf('\n*** Matrix %s \n', filename);
%     file = fullfile(matdir, [filename '.mat']);
%     benchmark_all_mat(mat_data, lambda, p, mat_rank, true, file);

    fprintf('\n*** Tensor %s \n', filename);
    file = fullfile(tendir, [filename '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add all other algorithms

%     % t-SVD
%     data.X = X; data.A = X;
%     params_tsvd.lambda = (0.1:0.05:0.5) ./ sqrt(max(size(X(:,:,1))));
%     params_tsvd.p = 1; params_tsvd.q = 1;
%     results.tsvd = benchmark(@geo_tensor_rpca, data, params_tsvd, true);
%     
%     % TNN
%     data.X = X; data.A = X;
%     params_tnn.lambda = (0.1:0.05:0.5) ./ sqrt(size(X,3)*max( size(X(:,:,1)) ) );
%     params_tnn.p = 1; params_tnn.q = 1;
%     results.tnn = benchmark(@geo_trpca_tnn, data, params_tnn, true);
%     
%     % 2D RPCA
%     data.X = X; data.A = X;
%     params_2d.lambda = 2*logspace(-1, -4, 20);
%     params_2d.p = 1; params_2d.q = 1;
%     
%     all_rpca2d = ...
%         {@rpca2d_l2, @rpca2d_l1, @rpca2d_gl_l2, @rpca2d_gl_l1};
%     all_rpca2d_names = ...
%         {'rpca2d_l2', 'rpca2d_l1', 'rpca2d_gl_l2', 'rpca2d_gl_l1'};
%     
%     for i=1:length(all_rpca2d)
%         nam_rpca2d = all_rpca2d_names{i};
%         fprintf('Benchmarking %s...\n', nam_rpca2d);
%         alg_rpca2d = @(X, opt)(geo_rpca2d_highway(all_rpca2d{i}, X, opt));
%         
%         results.(nam_rpca2d) = benchmark(alg_rpca2d, data, ...
%         params_2d, true);
%     end

    data.X = X; data.A = X;
    
    % Goldfarb
    lambda_opt = 1 / sqrt(max(size(X)));
    params_horpca.lambda = logspace(-1, 1, 30) * lambda_opt;
    params_horpca.p = 1; params_horpca.q = 1;
    results_part2.horpca_s = benchmark(@geo_horpca, data, params_horpca, true);
    
    % BRTF
    params_brtf.lambda = logspace(-7, 1, 30);
    params_brtf.p = 1; params_brtf.q = 1;
    results_part2.brtf = benchmark(@geo_brtf, data, params_brtf, true);
    
    save(file, 'results_part2', '-v7.3');
end

% save parameters
paramsfile = fullfile(rootdir, 'params_no_noise_part2.mat');
save(paramsfile, 'datafile', 'numims', 'lambda', 'p',...
    'mat_rank', 'ten_rank', 'ten_nrank_perc', ...
    'params_horpca', 'params_brtf');
