% quick_benchmark_mat_missing.m
% A quick benchmark on missing data for all matrix methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clc;
clear;
close all;

%% load data
datapath = '../data/';
%datafile = 'yaleb10.mat';
%datafile = 'highway.mat';
datafile = 'multipie10.mat';
load([datapath datafile]);

%% preprocess data
numims = [];
missing_perc = 0.1;
options.verbose = true;
options.return_error = false;

range = [min(X(:)), max(X(:))];
sizeX = size(X);
imsize = sizeX(rowidx);

X = arr2mat(X, rowidx, colidx);
if ~isempty(numims)
    X = X(:, 1:numims);
end

missing_ims = rand(1, size(X,2)) < missing_perc;
known_mask = true(size(X));
known_mask(:, missing_ims) = false;
options.error_weights = known_mask;
Y = X;
Y(~known_mask) = 0;

%% RPCA
options.lambda = 10;
options.method = 'alm';
tic
[Yr, info_r] = rpca(Y, options);
time_r = toc;
err_r = norm(X - Yr.A, 'fro') / norm(X, 'fro');
fprintf('RPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_r);
fprintf('\t error = %g \n', err_r);
fprintf('\t iterations = %d \n\n', info_r.iter(end));

disp_imdata(cat(3, X, Y, Yr.A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% BRPCA
options.lambda = 10;
options.numComp = 200;
tic
[Yb, info_b] = brpca(Y, options);
time_b = toc;
err_b = norm(X - Yb.A, 'fro') / norm(X, 'fro');
fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_b);
fprintf('\t error = %g \n', err_b);
fprintf('\t iterations = %d \n\n', info_b.iter(end));

disp_imdata(cat(3, X, Y, Yb.A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% IRPCA
options.lambda = 10;
options.method = 'lin';
tic
[Yi, info_i] = irpca(Y, options);
time_i = toc;
err_i = norm(X - Yi.A, 'fro') / norm(X, 'fro');
fprintf('IRPCA (%s): \n', options.method);
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_i);
fprintf('\t error = %g \n', err_i);
fprintf('\t iterations = %d \n\n', info_i.iter(end));

disp_imdata(cat(3, X, Y, Yi.A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% ORPCA
options.lambda = 10;
options.numComp = 100;
tic
[Yo, info_o] = orpca(Y, options);
time_o = toc;
err_o = norm(X - Yo.A, 'fro') / norm(X, 'fro');
fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_o);
fprintf('\t error = %g \n', err_o);
fprintf('\t iterations = %d \n\n', info_o.iter(end));

disp_imdata(cat(3, X, Y, Yo.A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

%% ROSL
options.lambda = 10;
options.numComp = 100;
%options.initMu = 1.0e-1;
tic
[Ys, info_s] = rosl(Y, options);
time_s = toc;
err_s = norm(X - Ys.A, 'fro') / norm(X, 'fro');
fprintf('ROSL: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_s);
fprintf('\t error = %g \n', err_s);
fprintf('\t iterations = %d \n\n', info_s.iter(end));

disp_imdata(cat(3, X, Y, Ys.A), imsize, {'Original', 'Corrupted', 'Low-rank'}, [1 3], range);

