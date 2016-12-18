% experiment_hall_no_noise.m
% Executes the background subtraction experiment on the hall video,
% without the addition of noise.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% directory where to write
rootdir = '/vol/bitbucket/mb2215/experiments/hall';
matdir = fullfile(rootdir, 'matrices');
tendir = fullfile(rootdir, 'tensors');
mkdir(matdir);
mkdir(tendir);

% data parameters
datapath = fullfile('..', 'data');
datafile = 'hall.mat';
numims = 300;
% load(fullfile(datapath, datafile));

% MODIFICATION BY MEHDI BAHRI
load_hall
X = O; clear O;

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

    filename = sprintf('%dx%dx%d_no_noise_rpca', size(thisX));

     fprintf('\n*** Matrix %s \n', filename);
     file = fullfile(matdir, [filename '.mat']);
%     benchmark_all_mat(mat_data, lambda, p, mat_rank, true, file);

%    fprintf('\n*** Tensor %s \n', filename);
%    file = fullfile(tendir, [filename '.mat']);
    % MODIFICATIONS BY MEHDI BAHRI
%     results = benchmark_all_ten(ten_data, lambda, p, ...
%         ten_nrank, ten_rank, true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add all other algorithms

    % RPCA
    data.X = X; data.A = X;
    params_rpca.lambda = 2*logspace(-1, -4, 20);
    params_rpca.p = 1; params_rpca.q = 1;
    results_rpca = benchmark(@geo_rpca, data, params_rpca, true);
    
    save(file, 'results_rpca', '-v7.3');
end

% save parameters
paramsfile = fullfile(rootdir, 'params_no_noise_rpca.mat');
save(paramsfile, 'datafile', 'numims', 'lambda', 'p',...
    'mat_rank', 'ten_rank', 'ten_nrank_perc', ...
    'params_rpca');
