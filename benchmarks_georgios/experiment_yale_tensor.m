% experiment_yale_tensor.m
% Executes the denoising experiment on Yale B, treating the data as a
% 4-mode tensor (with ID being a mode).
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% directory where to write
rootdir = '/vol/bitbucket/gp1813/experiments/yale';
tendir = fullfile(rootdir, 'tensors');
datadir = fullfile(rootdir, 'data');
mkdir(tendir);
mkdir(datadir);

% data parameters
datapath = fullfile('..', 'data');
datafile = 'yaleb10_tensor.mat';
noise_perc = [0.1 0.3 0.6 0.7];
patch_size = [10 30 40];
load(fullfile(datapath, datafile));

% algorithm parameters
lambda = logspace(-1, -4, 7);
p = 1;
ten_rank = 10 : 10 : 50;
ten_nrank_perc = 0.2 : 0.2 : 1.0;
ten_nrank = cell(1, length(ten_nrank_perc));
for i = 1:length(ten_nrank_perc)
    ten_nrank{i} = floor(ten_nrank_perc(i) * size(X));
end

ten_data.A = X;

% no noise
ten_data.X = X;

filename = sprintf('%dx%dx%dx%d_no_noise', size(X));

fprintf('\n*** Tensor %s \n', filename);
file = fullfile(tendir, [filename '.mat']);
benchmark_all_ten(ten_data, lambda, p, ten_nrank, ten_rank, true, file);

% salt and pepper noise
for np = 1:length(noise_perc)

    thisX_noise = add_noise(X, 'salt_pepper', noise_perc(np), rowidx);
    ten_data.X = thisX_noise;

    filename = sprintf('%dx%dx%dx%d_salt_pepper_%g', size(X), noise_perc(np));

    % save noisy data
    fprintf('\n*** Data %s \n', filename);
    file = fullfile(datadir, [filename '.mat']);
    save(file, 'thisX_noise');

    % run tensor methods
    fprintf('\n*** Tensor %s \n', filename);
    file = fullfile(tendir, [filename '.mat']);
    benchmark_all_ten(ten_data, lambda, p, ten_nrank, ten_rank, true, file);

end

% random patch
for ps = 1:length(patch_size)

    thisX_noise = add_noise(X, 'patch', patch_size(ps), rowidx);
    ten_data.X = thisX_noise;

    filename = sprintf('%dx%dx%dx%d_random_patch_%g', size(X), patch_size(ps));

    % save noisy data
    fprintf('\n*** Data %s \n', filename);
    file = fullfile(datadir, [filename '.mat']);
    save(file, 'thisX_noise');

    % run tensor methods
    fprintf('\n*** Tensor %s \n', filename);
    file = fullfile(tendir, [filename '.mat']);
    benchmark_all_ten(ten_data, lambda, p, ten_nrank, ten_rank, true, file);

end

% save parameters
paramsfile = fullfile(rootdir, 'params_tensor.mat');
save(paramsfile, 'datafile', 'noise_perc', 'patch_size', 'lambda', 'p', 'ten_rank', 'ten_nrank_perc');
