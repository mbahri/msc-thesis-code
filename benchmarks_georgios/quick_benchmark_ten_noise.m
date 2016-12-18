% quick_benchmark_ten_noise.m
% A quick benchmark on denoising for all tensor methods.
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
noise_type = 'salt_pepper';
noise_param = 0.10;
options.verbose = true;

Y = add_noise(X, noise_type, noise_param, rowidx);
range = [min(Y(:)), max(Y(:))];

sizeX = size(X);
imsize = sizeX(rowidx);

matX = arr2mat(X, rowidx, colidx);
matY = arr2mat(Y, rowidx, colidx);

%% RPCA
options.lambda = 0.005;
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
E = arr2mat(Yr.E, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A, E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% BRPCA
options.lambda = 0.01;
options.numComp = floor(0.8 * size(Y));
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
E = arr2mat(Yb.E, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A, E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% IRPCA
options.lambda = 0.001;
%options.numComp = [1 1 1];
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
E = arr2mat(Yi.E, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A, E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% ORPCA
options.lambda = 0.05;
options.numComp = floor(0.7 * size(Y));
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
E = arr2mat(Yo.E, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A, E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);

%% RCPD
options.lambda = 0.05;
options.method = 'sub';
options.numComp = 50;
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
E = arr2mat(Ycp.E, rowidx, colidx);
disp_imdata(cat(3, matX, matY, A, E), imsize, {'Original', 'Corrupted', 'Low-rank', 'Sparse'}, [1 4], range);
