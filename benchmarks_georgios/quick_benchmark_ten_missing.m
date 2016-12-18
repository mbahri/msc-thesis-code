% quick_benchmark_ten_missing.m
% A quick benchmark on missing data for all tensor methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clc;
clear;
close all;

%% create data
datapath = '../data/';
%datafile = 'yaleb10.mat';
%datafile = 'highway.mat';
datafile = 'multipie10.mat';
load([datapath datafile]);

%% preprocess data
missing_perc = 0.1;
options.verbose = true;
options.return_error = false;

range = [min(X(:)), max(X(:))];
sizeX = size(X);
imsize = sizeX(rowidx);

X = arr2mat(X, rowidx, colidx);
missing_ims = rand(1, size(X,2)) < missing_perc;
known_mask = true(size(X));
known_mask(:, missing_ims) = false;
Y = X;
Y(~known_mask) = 0;
X = mat2arr(X, sizeX, rowidx, colidx);
Y = mat2arr(Y, sizeX, rowidx, colidx);
options.error_weights = mat2arr(double(known_mask), sizeX, rowidx, colidx);

matX = arr2mat(X, rowidx, colidx);
matY = arr2mat(Y, rowidx, colidx);

%% RPCA
options.alpha = [1 1 1 1 0 0]; 
options.alpha = options.alpha / sum(options.alpha);
options.lambda = 10;
tic
[Yr, info_r] = ten_rpca(Y, options);
time_r = toc;
err_r = ew_norm(X - Yr.A, 2) / ew_norm(X, 2);
fprintf('RPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_r);
fprintf('\t error = %g \n', err_r);
fprintf('\t iterations = %d \n\n', info_r.iter(end));

A = arr2mat(Yr.A, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% BRPCA
options.alpha = [1 1 1 1 0 0]; 
options.alpha = options.alpha / sum(options.alpha);
options.lambda = 10;
options.numComp = floor([0.9 0.9 0.9 0.9 1 1] .* size(Y));
tic
[Yb, info_b] = ten_brpca(Y, options);
time_b = toc;
err_b = ew_norm(X - Yb.A, 2) / ew_norm(X, 2);
fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_b);
fprintf('\t error = %g \n', err_b);
fprintf('\t iterations = %d \n\n', info_b.iter(end));

A = arr2mat(Yb.A, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% IRPCA
options.alpha = [1 1 1 1 0 0]; 
options.alpha = options.alpha / sum(options.alpha);
options.lambda = 10;
options.numComp = floor(0.5 * size(Y));
options.method = 'sub';
tic
[Yi, info_i] = ten_irpca(Y, options);
time_i = toc;
err_i = ew_norm(X - Yi.A, 2) / ew_norm(X, 2);
fprintf('IRPCA (%s): \n', options.method);
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_i);
fprintf('\t error = %g \n', err_i);
fprintf('\t iterations = %d \n\n', info_i.iter(end));

A = arr2mat(Yi.A, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% ORPCA
options.alpha = [1 1 1 1 0 0]; 
options.alpha = options.alpha / sum(options.alpha);
options.lambda = 10;
options.numComp = floor([1 1 1 1 1 1] .* size(Y));
tic;
[Yo, info_o] = ten_orpca(Y, options);
time_o = toc;
err_o = ew_norm(X - Yo.A, 2) / ew_norm(X, 2);
fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_o);
fprintf('\t error = %g \n', err_o);
fprintf('\t iterations = %d \n\n', info_o.iter(end));

A = arr2mat(Yo.A, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% RCPD
options.alpha = [1 1 1 1 0 0]; 
options.alpha = options.alpha / sum(options.alpha);
options.lambda = 10;
options.method = 'sub';
options.numComp = 80;
tic;
[Ycp, info_cp] = ten_rcpd(Y, options);
time_cp = toc;
err_cp = ew_norm(X - Ycp.A, 2) / ew_norm(X, 2);
fprintf('RCPD (%s): \n', options.method);
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_cp);
fprintf('\t error = %g \n', err_cp);
fprintf('\t iterations = %d \n\n', info_cp.iter(end));

A = arr2mat(Ycp.A, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);
