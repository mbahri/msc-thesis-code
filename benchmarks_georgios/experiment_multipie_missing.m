% experiment_multipie_missing.m
% Executes the image reconstruction experiment on Multi-PIE.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% directory where to write
rootdir = '/vol/bitbucket/gp1813/experiments/multipie';
tendir = fullfile(rootdir, 'tensors');
mkdir(tendir);

% data parameters
datapath = fullfile('..', 'data');
datafile = 'multipie10.mat';
load(fullfile(datapath, datafile));
missing_perc = 0.2;

% algorithm parameters
alpha = [1 1 1 1 0 0];
lambda = [1 10 100 1000];
p = 1;
ten_rank = 40;
ten_nrank_perc = 1.0;
ten_nrank = cell(1, length(ten_nrank_perc));
for i = 1:length(ten_nrank_perc)
    ten_nrank{i} = floor(ten_nrank_perc(i) * size(X));
end

ten_data.A = X;

% remove images to be missing 
matX = arr2mat(X, rowidx, colidx);
missing_ims = rand(1, size(matX,2)) < missing_perc;
known_mask = true(size(matX));
known_mask(:, missing_ims) = false;
matX(~known_mask) = 0;
ten_data.X = mat2arr(matX, size(X), rowidx, colidx);
error_weights = mat2arr(double(known_mask), size(X), rowidx, colidx);

filename = sprintf('%dx%dx%dx%dx%dx%d_missing_%g', size(X), missing_perc);

fprintf('\n*** Tensor %s \n', filename);
file = fullfile(tendir, [filename '.mat']);
benchmark_all_ten(ten_data, lambda, p, ten_nrank, ten_rank, true, file, alpha, error_weights);

% save parameters
paramsfile = fullfile(rootdir, sprintf('params_missing_%g.mat', missing_perc));
save(paramsfile, 'datafile', 'missing_perc', 'alpha', 'lambda', 'p', 'ten_rank', 'ten_nrank_perc', 'error_weights');
