% experiment_yale_random_patch.m
% Executes the denoising experiment on Yale B, with addition of a random 
% patch.
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
patch_size = [10 30 40];
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
    
    % random patch
    for ps = 1:length(patch_size)
        
        thisX_noise = add_noise(thisX, 'patch', patch_size(ps), rowidx);
        mat_data.X = arr2mat(thisX_noise, rowidx, colidx);
        ten_data.X = thisX_noise;

        filename = sprintf('%dx%dx%d_random_patch_%g', size(thisX), patch_size(ps));
        
        % save noisy data
        fprintf('\n*** Data %s \n', filename);
        file = fullfile(datadir, [filename '.mat']);
        save(file, 'thisX_noise');

        % run matrix methods
        fprintf('\n*** Matrix %s \n', filename);
        file = fullfile(matdir, [filename '.mat']);
        benchmark_all_mat(mat_data, lambda, p, mat_rank, true, file);

        % run tensor methods
        fprintf('\n*** Tensor %s \n', filename);
        file = fullfile(tendir, [filename '.mat']);
        benchmark_all_ten(ten_data, lambda, p, ten_nrank, ten_rank, true, file);
    
    end
end

% save parameters
paramsfile = fullfile(rootdir, 'params_random_patch.mat');
save(paramsfile, 'datafile', 'numims', 'patch_size', 'lambda', 'p', 'mat_rank', 'ten_rank', 'ten_nrank_perc');
