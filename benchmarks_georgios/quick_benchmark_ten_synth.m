% quick_benchmark_ten_synth.m
% A quick benchmark on synthetic data for all tensor methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clc;
clear;
close all;

%% create data
opt.type = 'tucker_tensor';
opt.dims = [50 50 50];
opt.rank_perc = 0.05;
opt.data_magn = 1;
opt.err_perc = 0.1;
opt.err_magn = 100;
opt.return_error = false;

data = create_synthetic_data(opt);

%% RPCA
options.lambda = 0.05;
tic
[Yr, info_r] = ten_rpca(data.X, options);
time_r = toc;
err_r = norm(data.A - Yr.A) / norm(data.A);
fprintf('RPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_r);
fprintf('\t error = %g \n', err_r);
fprintf('\t iterations = %d \n\n', info_r.iter(end));

%% BRPCA
options.lambda = 0.05;
options.numComp = floor(opt.rank_perc * opt.dims);
tic
[Yb, info_b] = ten_brpca(data.X, options);
time_b = toc;
err_b = ew_norm(data.A - Yb.A, 2) / ew_norm(data.A, 2);
fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_b);
fprintf('\t error = %g \n', err_b);
fprintf('\t iterations = %d \n\n', info_b.iter(end));

%% IRPCA
options.lambda = 0.0001;
options.numComp = floor(opt.rank_perc * opt.dims);
options.method = 'sub';
tic
[Yi, info_i] = ten_irpca(data.X, options);
time_i = toc;
err_i = norm(data.A - Yi.A) / norm(data.A);
fprintf('IRPCA (%s): \n', options.method);
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_i);
fprintf('\t error = %g \n', err_i);
fprintf('\t iterations = %d \n\n', info_i.iter(end));

%% ORPCA
options.lambda = 0.05;
options.numComp = floor(opt.rank_perc * opt.dims);
tic;
[Yo, info_o] = ten_orpca(data.X, options);
time_o = toc;
err_o = norm(data.A - Yo.A) / norm(data.A);
fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_o);
fprintf('\t error = %g \n', err_o);
fprintf('\t iterations = %d \n\n', info_o.iter(end));

%% RCPD
options.lambda = 0.05;
options.method = 'sub';
options.numComp = floor(opt.rank_perc * min(opt.dims));
tic;
[Ycp, info_cp] = ten_rcpd(data.X, options);
time_cp = toc;
err_cp = norm(data.A - Ycp.A) / norm(data.A);
fprintf('RCPD (%s): \n', options.method);
fprintf('\t lambda = %g \n', options.lambda);
fprintf('\t time = %g sec \n', time_cp);
fprintf('\t error = %g \n', err_cp);
fprintf('\t iterations = %d \n\n', info_cp.iter(end));
