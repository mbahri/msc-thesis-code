% quick_benchmark_mat_noise.m
% A quick benchmark on denoising for all matrix methods.
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
noise_type = 'whole_pics';
noise_param = 0.1;
options.verbose = true;

Y = add_noise(X, noise_type, noise_param, rowidx);
range = [min(Y(:)), max(Y(:))];

sizeX = size(X);
imsize = sizeX(rowidx);

X = arr2mat(X, rowidx, colidx);
Y = arr2mat(Y, rowidx, colidx);
if ~isempty(numims)
    X = X(:, 1:numims);
    Y = Y(:, 1:numims);
end

%% CPCA
options.numComp = 20;
tic
[Yc, info_c] = cpca(Y, options);
time_c = toc;
err_c = norm(X - Yc.A, 'fro') / norm(X, 'fro');
fprintf('CPCA: \n');
fprintf('\t time = %g sec \n', time_c);
fprintf('\t error = %g \n', err_c);
fprintf('\t iterations = %d \n\n', info_c.iter(end));

disp_imdata(cat(3, X, Y, Yc.A, Yc.E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% RPCA
options.lambda = 0.02;
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

disp_imdata(cat(3, X, Y, Yr.A, Yr.E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% BRPCA
options.lambda = 0.05;
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

disp_imdata(cat(3, X, Y, Yb.A, Yb.E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% IRPCA
options.lambda = 0.05;
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

disp_imdata(cat(3, X, Y, Yi.A, Yi.E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% ORPCA
options.lambda = 0.05;
options.numComp = 50;
tic
[Yo, info_o] = orpca(Y, options);
time_o = toc;
err_o = norm(X - Yo.A, 'fro') / norm(X, 'fro');
fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_o);
fprintf('\t error = %g \n', err_o);
fprintf('\t iterations = %d \n\n', info_o.iter(end));

disp_imdata(cat(3, X, Y, Yo.A, Yo.E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% ROSL
options.lambda = 0.05;
options.numComp = 50;
options.initMu = 1.0e-1;
tic
[Ys, info_s] = rosl(Y, options);
time_s = toc;
err_s = norm(X - Ys.A, 'fro') / norm(X, 'fro');
fprintf('ROSL: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_s);
fprintf('\t error = %g \n', err_s);
fprintf('\t iterations = %d \n\n', info_s.iter(end));

disp_imdata(cat(3, X, Y, Ys.A, Ys.E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

