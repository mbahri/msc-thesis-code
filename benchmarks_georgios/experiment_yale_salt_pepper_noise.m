% experiment_yale_salt_pepper_noise.m
% Executes the denoising experiment on Yale B, with addition of salt &
% pepper noise.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% directory where to write
rootdir = '/vol/bitbucket/mb2215/experiments/yale';
matdir = fullfile(rootdir, 'matrices');
tendir = fullfile(rootdir, 'tensors');
datadir = fullfile(rootdir, 'data');
mkdir(matdir);
mkdir(tendir);
mkdir(datadir);

% data parameters
datapath = fullfile('..', 'data');
datafile = 'yaleb10.mat';
numims = [64 640];
noise_perc = [0.1 0.3 0.6 0.7];

% load(fullfile(datapath, datafile));
load_yale;
X = O;
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
        
    % salt and pepper noise
    for np = 1:length(noise_perc)
        
        thisX_noise = add_noise(thisX, 'salt_pepper', noise_perc(np), rowidx);
        mat_data.X = arr2mat(thisX_noise, rowidx, colidx);
        ten_data.X = thisX_noise;

        filename = sprintf('%dx%dx%d_salt_pepper_%g', size(thisX), noise_perc(np));

        % save noisy data
        fprintf('\n*** Data %s \n', filename);
        file = fullfile(datadir, [filename '.mat']);
        save(file, 'thisX_noise');

%         % run matrix methods
%         fprintf('\n*** Matrix %s \n', filename);
%         file = fullfile(matdir, [filename '.mat']);
%         benchmark_all_mat(mat_data, lambda, p, mat_rank, true, file);

        % run tensor methods
        fprintf('\n*** Tensor %s \n', filename);
        file = fullfile(tendir, [filename '.mat']);
        results = benchmark_all_ten(ten_data, lambda, p, ten_nrank, ...
            ten_rank, true);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add all other algorithms

        % t-SVD
        data.X = X; data.A = X;
        params_tsvd.lambda = (0.5:0.05:1.5) ./ sqrt(max( size(X(:,:,1)) ));
        params_tsvd.p = 1; params_tsvd.q = 1;
        results.tsvd = benchmark(@geo_tensor_rpca, data, params_tsvd, true);

        % TNN
        data.X = X; data.A = X;
        params_tnn.lambda = (0.5:0.05:1.5) ./ sqrt(size(X,3)*max( size(X(:,:,1)) ));
        params_tnn.p = 1; params_tnn.q = 1;
        results.tnn = benchmark(@geo_trpca_tnn, data, params_tnn, true);

        % 2D RPCA
        data.X = X; data.A = X;
        params_2d.lambda = (0.5:5e-2:1.5) ./ sqrt(size(X,3)*max( size(X(:,:,1)) ));
        params_2d.p = 1; params_2d.q = 1;
        results.rpca2d_l2 = benchmark(@geo_rpca2d_highway, data, ...
            params_2d, true);
        results.rpca2d_l1 = benchmark(@geo_rpca2d_highway, data, ...
            params_2d, true);

        % Goldfarb
        % TODO

        save(file, 'results');

    end
end

% save parameters
paramsfile = fullfile(rootdir, 'params_salt_pepper.mat');
save(paramsfile, 'datafile', 'numims', 'noise_perc', 'lambda', 'p', ...
    'mat_rank', 'ten_rank', 'ten_nrank_perc',...
    'params_2d', 'params_tsvd', 'params_tnn');
