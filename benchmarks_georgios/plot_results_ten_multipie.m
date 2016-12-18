% plot_results_ten_multipie.m
% Plots the results of the image reconstruction experiment on Multi-PIE
% using tensor methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
close all;
clc;

% choose data
dims = [10 5 6 5 40 30];
missing_perc = 0.2;
method_type = 'tensors';

% load data
path = '/vol/bitbucket/gp1813/experiments/multipie';
file = sprintf('%dx%dx%dx%dx%dx%d_missing_%g.mat', dims, missing_perc);
load(fullfile(path, method_type, file));
load(fullfile(path, sprintf('params_missing_%g.mat', missing_perc)));

% plot parameters
axis_fontsize = 14;
linewidth = 3;
markersize_big = 10;
markersize_sml = 5;
colour = distinguishable_colors(7);
colour = colour([1 5 4 3 7 6 2], :);
marker = {'d', 'x', 's', '>', '<', '+', '*'};
linestyle = ':';
position = [294 582 1051 392];

% display pics params
load('../data/multipie10.mat');
imsize = dims(rowidx);
range = [min(X(:)) max(X(:))];

% figure save path
savepath = fullfile('..', 'results', 'multipie');

% squeeze matrices
results.brpca.err  = squeeze(results.brpca.err );
results.brpca.time = squeeze(results.brpca.time);
results.brpca.info = squeeze(results.brpca.info);
results.brpca.A    = squeeze(results.brpca.A   );

results.irpca.sub.err  = squeeze(results.irpca.sub.err );
results.irpca.sub.time = squeeze(results.irpca.sub.time);
results.irpca.sub.info = squeeze(results.irpca.sub.info);
results.irpca.sub.A    = squeeze(results.irpca.sub.A   );

results.irpca.lin.err  = squeeze(results.irpca.lin.err );
results.irpca.lin.time = squeeze(results.irpca.lin.time);
results.irpca.lin.info = squeeze(results.irpca.lin.info);
results.irpca.lin.A    = squeeze(results.irpca.lin.A   );

results.orpca.err  = squeeze(results.orpca.err );
results.orpca.time = squeeze(results.orpca.time);
results.orpca.info = squeeze(results.orpca.info);
results.orpca.A    = squeeze(results.orpca.A   );

results.rcpd.sub.err  = squeeze(results.rcpd.sub.err );
results.rcpd.sub.time = squeeze(results.rcpd.sub.time);
results.rcpd.sub.info = squeeze(results.rcpd.sub.info);
results.rcpd.sub.A    = squeeze(results.rcpd.sub.A   );

results.rcpd.lin.err  = squeeze(results.rcpd.lin.err );
results.rcpd.lin.time = squeeze(results.rcpd.lin.time);
results.rcpd.lin.info = squeeze(results.rcpd.lin.info);
results.rcpd.lin.A    = squeeze(results.rcpd.lin.A   );

%% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_corr  = containers.Map;
best_A     = containers.Map;

missingX = select_best_missing({X}, X, error_weights, rowidx, colidx);

[best_A('rpca'     ), best_corr('rpca'     ), best_l_idx('rpca'     )                         ] = select_best_missing(results. rpca.    A, X, error_weights, rowidx, colidx);
[best_A('brpca'    ), best_corr('brpca'    ), best_l_idx('brpca'    ), best_r_idx('brpca'    )] = select_best_missing(results.brpca.    A, X, error_weights, rowidx, colidx);
[best_A('irpca_sub'), best_corr('irpca_sub'), best_l_idx('irpca_sub'), best_r_idx('irpca_sub')] = select_best_missing(results.irpca.sub.A, X, error_weights, rowidx, colidx);
[best_A('irpca_lin'), best_corr('irpca_lin'), best_l_idx('irpca_lin'), best_r_idx('irpca_lin')] = select_best_missing(results.irpca.lin.A, X, error_weights, rowidx, colidx);
[best_A('orpca'    ), best_corr('orpca'    ), best_l_idx('orpca'    ), best_r_idx('orpca'    )] = select_best_missing(results.orpca.    A, X, error_weights, rowidx, colidx);
[best_A('rcpd_sub' ), best_corr('rcpd_sub' ), best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )] = select_best_missing(results. rcpd.sub.A, X, error_weights, rowidx, colidx);
[best_A('rcpd_lin' ), best_corr('rcpd_lin' ), best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )] = select_best_missing(results. rcpd.lin.A, X, error_weights, rowidx, colidx);

maxiter = max([results. rpca.    info{best_l_idx('rpca'     )                         }.iter(end) ...
               results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca'    )}.iter(end) ...
               results.irpca.sub.info{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}.iter(end) ...
               results.irpca.lin.info{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}.iter(end) ...
               results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca'    )}.iter(end) ...
               results. rcpd.sub.info{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}.iter(end) ...
               results. rcpd.lin.info{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}.iter(end) ] );

%% time plots
cnt = 1;

fig = figure('Position', position);
semilogx(lambda, results. rpca.    time,                            linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogx(lambda, results.brpca.    time(:,best_r_idx('brpca'    )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.sub.time(:,best_r_idx('irpca_sub')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.lin.time(:,best_r_idx('irpca_lin')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.orpca.    time(:,best_r_idx('orpca'    )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results. rcpd.sub.time(:,best_r_idx('rcpd_sub' )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results. rcpd.lin.time(:,best_r_idx('rcpd_lin' )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Total execution time [sec]', 'FontSize', axis_fontsize);
xlabel('\lambda', 'FontSize', axis_fontsize);
legend('RPCA', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'RHOSVD', 'RCPD (sub)', 'RCPD (lin)', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([min(lambda) max(lambda)]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_time.eps'), 'psc2');

%% performance plots
cnt = 1;

fig = figure('Position', position);
semilogx(lambda, results. rpca.    err,                            linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogx(lambda, results.brpca.    err(:,best_r_idx('brpca'    )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.sub.err(:,best_r_idx('irpca_sub')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.lin.err(:,best_r_idx('irpca_lin')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.orpca.    err(:,best_r_idx('orpca'    )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results. rcpd.sub.err(:,best_r_idx('rcpd_sub' )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results. rcpd.lin.err(:,best_r_idx('rcpd_lin' )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Relative recovery error',   'FontSize', axis_fontsize);
xlabel('\lambda', 'FontSize', axis_fontsize);
legend('RPCA', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'RHOSVD', 'RCPD (sub)', 'RCPD (lin)', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([min(lambda) max(lambda)]);
%ylim([1.0e-7 1.0e+2]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_performance.eps'), 'psc2');

%% convergence plots
iter_show = maxiter; cnt = 1;

fig = figure('Position', position); 
semilogy(results. rpca.    info{best_l_idx('rpca'     )                         }.iter, results. rpca.    info{best_l_idx('rpca'     )                         }.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogy(results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca'    )}.iter, results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca')    }.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.sub.info{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}.iter, results.irpca.sub.info{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.lin.info{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}.iter, results.irpca.lin.info{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca'    )}.iter, results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca'    )}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results. rcpd.sub.info{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}.iter, results. rcpd.sub.info{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results. rcpd.lin.info{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}.iter, results. rcpd.lin.info{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Convergence criterion', 'FontSize', axis_fontsize);
xlabel('Number of iterations', 'FontSize', axis_fontsize);
legend('RPCA', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'RHOSVD', 'RCPD (sub)', 'RCPD (lin)', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([1 iter_show]);
ylim([1.0e-7, 1]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_convergence.eps'), 'psc2');


%% best results
fprintf('RPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca')));
fprintf('\t time = %g sec \n', results.rpca.time(best_l_idx('rpca')));
fprintf('\t correlation = %g \n', best_corr('rpca'));
fprintf('\t iterations = %d \n\n', results.rpca.info{best_l_idx('rpca')}.iter(end));

fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('brpca')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('brpca')));
fprintf('\t time = %g sec \n', results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')));
fprintf('\t correlation = %g \n', best_corr('brpca'));
fprintf('\t iterations = %d \n\n', results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end));

fprintf('IRPCA (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_sub')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('irpca_sub')));
fprintf('\t time = %g sec \n', results.irpca.sub.time(best_l_idx('irpca_sub'),best_r_idx('irpca_sub')));
fprintf('\t correlation = %g \n', best_corr('irpca_sub'));
fprintf('\t iterations = %d \n\n', results.irpca.sub.info{best_l_idx('irpca_sub'),best_r_idx('irpca_sub')}.iter(end));

fprintf('IRPCA (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_lin')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('irpca_lin')));
fprintf('\t time = %g sec \n', results.irpca.lin.time(best_l_idx('irpca_lin'),best_r_idx('irpca_lin')));
fprintf('\t correlation = %g \n', best_corr('irpca_lin'));
fprintf('\t iterations = %d \n\n', results.irpca.lin.info{best_l_idx('irpca_lin'),best_r_idx('irpca_lin')}.iter(end));

fprintf('RHOSVD: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('orpca')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('orpca')));
fprintf('\t time = %g sec \n', results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')));
fprintf('\t correlation = %g \n', best_corr('orpca'));
fprintf('\t iterations = %d \n\n', results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end));

fprintf('RCPD (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rcpd_sub')));
fprintf('\t rank = %g \n', ten_rank(best_r_idx('rcpd_sub')));
fprintf('\t time = %g sec \n', results.rcpd.sub.time(best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')));
fprintf('\t correlation = %g \n', best_corr('rcpd_sub'));
fprintf('\t iterations = %d \n\n', results.rcpd.sub.info{best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')}.iter(end));

fprintf('RCPD (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rcpd_lin')));
fprintf('\t rank = %g \n', ten_rank(best_r_idx('rcpd_lin')));
fprintf('\t time = %g sec \n', results.rcpd.lin.time(best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')));
fprintf('\t correlation = %g \n', best_corr('rcpd_lin'));
fprintf('\t iterations = %d \n\n', results.rcpd.lin.info{best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')}.iter(end));

%% print as latex table
fprintf('RPCA & $%.3f$ & $%.3f$ & $%d$ & $%g$ \\\\ \n', ...
    best_corr('rpca'), ...
    results.rpca.time(best_l_idx('rpca')), ...
    results.rpca.info{best_l_idx('rpca')}.iter(end), ...
    lambda(best_l_idx('rpca')) );

fprintf('BRPCA & $%.3f$ & $%.3f$ & $%d$ & $%g$  \\\\ \n', ...
    best_corr('brpca'), ...
    results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')), ...
    results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end), ...
    lambda(best_l_idx('brpca')) );

fprintf('IRPCA (sub) & $%.3f$ & $%.3f$ & $%d$ & $%g$ \\\\ \n', ...
    best_corr('irpca_sub'), ...
    results.irpca.sub.time(best_l_idx('irpca_sub'),best_r_idx('irpca_sub')), ...
    results.irpca.sub.info{best_l_idx('irpca_sub'),best_r_idx('irpca_sub')}.iter(end), ...
    lambda(best_l_idx('irpca_sub')) );

fprintf('IRPCA (lin) & $%.3f$ & $%.3f$ & $%d$ & $%g$ \\\\ \n', ...
    best_corr('irpca_lin'), ...
    results.irpca.lin.time(best_l_idx('irpca_lin'),best_r_idx('irpca_lin')), ...
    results.irpca.lin.info{best_l_idx('irpca_lin'),best_r_idx('irpca_lin')}.iter(end), ...
    lambda(best_l_idx('irpca_lin')) );

fprintf('RHOSVD & $%.3f$ & $%.3f$ & $%d$ & $%g$ \\\\ \n', ...
    best_corr('orpca'), ...
    results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')), ...
    results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end), ...
    lambda(best_l_idx('orpca')) );

fprintf('RCPD (sub) & $%.3f$ & $%.3f$ & $%d$ & $%g$ \\\\ \n', ...
    best_corr('rcpd_sub'), ...
    results.rcpd.sub.time(best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')), ...
    results.rcpd.sub.info{best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')}.iter(end), ...
    lambda(best_l_idx('rcpd_sub')) );

fprintf('RCPD (lin) & $%.3f$ & $%.3f$ & $%d$ & $%g$ \\\\ \n', ...
    best_corr('rcpd_lin'), ...
    results.rcpd.lin.time(best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')), ...
    results.rcpd.lin.info{best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')}.iter(end), ...
    lambda(best_l_idx('rcpd_lin')) );

%% display images
disp_imdata(cat(3, missingX, ...
    best_A('rpca'     ), ...
    best_A('brpca'    ), ...
    best_A('irpca_sub'), ...
    best_A('irpca_lin'), ...
    best_A('orpca'    ), ...
    best_A('rcpd_sub' ), ...
    best_A('rcpd_lin' ) ), ...
    imsize, ...
    {'Original', 'RPCA', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'RHOSVD', 'RCPD (sub)', 'RCPD (lin)'}, ...
    [2 4], []);

%% save images
switch missing_perc
    case 0.1
        im_idx = [113 33 59 11 26];
    case 0.2
        im_idx = [258 72 141 25 61];
    case 0.5
        im_idx = [609 170 341 55 150];
    case 0.8
        im_idx = [992 265 552 93 234];
    otherwise
        error('Image indices not specified yet.');
end
imformat = 'png';

for i = 1:length(im_idx)
    imfile = sprintf('ten_multipie_%g_original_%d.%s' , missing_perc, im_idx(i), imformat); im = missingX           ; im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_rpca_%d.%s'     , missing_perc, im_idx(i), imformat); im = best_A('rpca'     ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_brpca_%d.%s'    , missing_perc, im_idx(i), imformat); im = best_A('brpca'    ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_irpca_sub_%d.%s', missing_perc, im_idx(i), imformat); im = best_A('irpca_sub'); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_irpca_lin_%d.%s', missing_perc, im_idx(i), imformat); im = best_A('irpca_lin'); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_orpca_%d.%s'    , missing_perc, im_idx(i), imformat); im = best_A('orpca'    ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_rcpd_sub_%d.%s' , missing_perc, im_idx(i), imformat); im = best_A('rcpd_sub' ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_multipie_%g_rcpd_lin_%d.%s' , missing_perc, im_idx(i), imformat); im = best_A('rcpd_lin' ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
end
