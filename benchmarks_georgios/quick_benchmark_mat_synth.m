% quick_benchmark_mat_synth.m
% A quick benchmark on synthetic data for all matrix methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clc;
clear;
close all;

%% create data
opt.type = 'matrix';
opt.dims = [500 500];
opt.rank_perc = 0.05;
opt.data_magn = 1;
opt.err_perc = 0.1;
opt.err_magn = 100;
opt.return_error = false;

data = create_synthetic_data(opt);

%% CPCA
options.numComp = floor(opt.rank_perc * min(opt.dims));
tic
[Yc, info_c] = cpca(data.X, options);
time_c = toc;
err_c = norm(data.A - Yc.A, 'fro') / norm(data.A, 'fro');
fprintf('CPCA: \n');
fprintf('\t time = %g sec \n', time_c);
fprintf('\t error = %g \n', err_c);
fprintf('\t iterations = %d \n\n', info_c.iter(end));

%% RPCA
options.lambda = 0.05;
options.method = 'alm';
tic
[Yr, info_r] = rpca(data.X, options);
time_r = toc;
err_r = norm(data.A - Yr.A, 'fro') / norm(data.A, 'fro');
fprintf('RPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_r);
fprintf('\t error = %g \n', err_r);
fprintf('\t iterations = %d \n\n', info_r.iter(end));

%% BRPCA
options.lambda = 0.05;
options.numComp = floor(opt.rank_perc * min(opt.dims));
tic
[Yb, info_b] = brpca(data.X, options);
time_b = toc;
err_b = norm(data.A - Yb.A, 'fro') / norm(data.A, 'fro');
fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_b);
fprintf('\t error = %g \n', err_b);
fprintf('\t iterations = %d \n\n', info_b.iter(end));

%% IRPCA
options.lambda = 0.05;
options.method = 'lin';
tic
[Yi, info_i] = irpca(data.X, options);
time_i = toc;
err_i = norm(data.A - Yi.A, 'fro') / norm(data.A, 'fro');
fprintf('IRPCA (%s): \n', options.method);
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_i);
fprintf('\t error = %g \n', err_i);
fprintf('\t iterations = %d \n\n', info_i.iter(end));

%% ORPCA
options.lambda = 0.05;
options.numComp = floor(opt.rank_perc * min(opt.dims));
tic
[Yo, info_o] = orpca(data.X, options);
time_o = toc;
err_o = norm(data.A - Yo.A, 'fro') / norm(data.A, 'fro');
fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_o);
fprintf('\t error = %g \n', err_o);
fprintf('\t iterations = %d \n\n', info_o.iter(end));

%% ROSL
options.lambda = 0.05;
options.numComp = floor(opt.rank_perc * min(opt.dims));
options.initMu = 1.0e-1;
tic
[Ys, info_s] = rosl(data.X, options);
time_s = toc;
err_s = norm(data.A - Ys.A, 'fro') / norm(data.A, 'fro');
fprintf('ROSL: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_s);
fprintf('\t error = %g \n', err_s);
fprintf('\t iterations = %d \n\n', info_s.iter(end));
