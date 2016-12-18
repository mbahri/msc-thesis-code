% plot_results_mat_highway.m
% Plots the results of the background subtraction experiment on the highway
% video using matrix methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
close all;
clc;

% choose data
dims = [48 64 400];
method_type = 'matrices';

% load data
path = '/vol/bitbucket/gp1813/experiments/highway';
file = sprintf('%dx%dx%d_no_noise.mat', dims);
load(fullfile(path, method_type, file));
load(fullfile(path, 'params_no_noise.mat'));

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
load('../data/highway.mat');
imsize = dims(rowidx);
range = [min(X(:)) max(X(:))];
X  = arr2mat(X,  rowidx, colidx);
GT = arr2mat(GT, rowidx, colidx);

% figure save path
savepath = fullfile('..', 'results', 'highway');

% squeeze matrices
results.brpca.err  = squeeze(results.brpca.err );
results.brpca.time = squeeze(results.brpca.time);
results.brpca.info = squeeze(results.brpca.info);
results.brpca.A    = squeeze(results.brpca.A   );

results.orpca.err  = squeeze(results.orpca.err );
results.orpca.time = squeeze(results.orpca.time);
results.orpca.info = squeeze(results.orpca.info);
results.orpca.A    = squeeze(results.orpca.A   );

results.rosl.err  = squeeze(results.rosl.err );
results.rosl.time = squeeze(results.rosl.time);
results.rosl.info = squeeze(results.rosl.info);
results.rosl.A    = squeeze(results.rosl.A   );

%% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_tpr   = containers.Map;
best_fpr   = containers.Map;
best_auc   = containers.Map;
best_foreground = containers.Map;

true_mask = GT > min(GT(:)) + 0.1;
threshold = linspace(0, 1, 101);

[best_foreground('rpca_alm' ), best_tpr('rpca_alm' ), best_fpr('rpca_alm' ), best_auc('rpca_alm' ), best_l_idx('rpca_alm' )                     ] = select_best_foreground(results. rpca.alm.A, X, true_mask, threshold);
[best_foreground('rpca_apg' ), best_tpr('rpca_apg' ), best_fpr('rpca_apg' ), best_auc('rpca_apg' ), best_l_idx('rpca_apg' )                     ] = select_best_foreground(results. rpca.apg.A, X, true_mask, threshold);
[best_foreground('brpca'    ), best_tpr('brpca'    ), best_fpr('brpca'    ), best_auc('brpca'    ), best_l_idx('brpca'    ), best_r_idx('brpca')] = select_best_foreground(results.brpca.    A, X, true_mask, threshold);
[best_foreground('irpca_sub'), best_tpr('irpca_sub'), best_fpr('irpca_sub'), best_auc('irpca_sub'), best_l_idx('irpca_sub')                     ] = select_best_foreground(results.irpca.sub.A, X, true_mask, threshold);
[best_foreground('irpca_lin'), best_tpr('irpca_lin'), best_fpr('irpca_lin'), best_auc('irpca_lin'), best_l_idx('irpca_lin')                     ] = select_best_foreground(results.irpca.lin.A, X, true_mask, threshold);
[best_foreground('orpca'    ), best_tpr('orpca'    ), best_fpr('orpca'    ), best_auc('orpca'    ), best_l_idx('orpca'    ), best_r_idx('orpca')] = select_best_foreground(results.orpca.    A, X, true_mask, threshold);
[best_foreground('rosl'     ), best_tpr('rosl'     ), best_fpr('rosl'     ), best_auc('rosl'     ), best_l_idx('rosl'     ), best_r_idx('rosl' )] = select_best_foreground(results. rosl.    A, X, true_mask, threshold);

maxiter = max([results. rpca.alm.info{best_l_idx('rpca_alm' )                     }.iter(end) ...
               results. rpca.apg.info{best_l_idx('rpca_apg' )                     }.iter(end) ...
               results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca')}.iter(end) ...
               results.irpca.sub.info{best_l_idx('irpca_sub')                     }.iter(end) ...
               results.irpca.lin.info{best_l_idx('irpca_lin')                     }.iter(end) ...
               results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca')}.iter(end) ...
               results. rosl.    info{best_l_idx('rosl'     ), best_r_idx('rosl' )}.iter(end)] );

%% time plots
cnt = 1;

fig = figure('Position', position);
loglog(lambda, results. rpca.alm.time,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
loglog(lambda, results. rpca.apg.time,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.brpca.    time(:,best_r_idx('brpca')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.irpca.sub.time,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.irpca.lin.time,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.orpca.    time(:,best_r_idx('orpca')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results. rosl.    time(:,best_r_idx('rosl' )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Total execution time [sec]', 'FontSize', axis_fontsize);
xlabel('\lambda', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([min(lambda) max(lambda)]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_time.eps'), 'psc2');

%% performance plots
cnt = 1;

fig = figure('Position', position);
semilogx(lambda, results. rpca.alm.err,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogx(lambda, results. rpca.apg.err,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.brpca.    err(:,best_r_idx('brpca')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.sub.err,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.lin.err,                        linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.orpca.    err(:,best_r_idx('orpca')), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results. rosl.    err(:,best_r_idx('rosl' )), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Relative recovery error',   'FontSize', axis_fontsize);
xlabel('\lambda', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([min(lambda) max(lambda)]);
%ylim([1.0e-7 1.0e+2]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_performance.eps'), 'psc2');

%% convergence plots
iter_show = maxiter; cnt = 1;

fig = figure('Position', position); 
semilogy(results. rpca.alm.info{best_l_idx('rpca_alm' )                     }.iter, results. rpca.alm.info{best_l_idx('rpca_alm' )                     }.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogy(results. rpca.apg.info{best_l_idx('rpca_apg' )                     }.iter, results. rpca.apg.info{best_l_idx('rpca_apg' )                     }.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca')}.iter, results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca')}.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.sub.info{best_l_idx('irpca_sub')                     }.iter, results.irpca.sub.info{best_l_idx('irpca_sub')                     }.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.lin.info{best_l_idx('irpca_lin')                     }.iter, results.irpca.lin.info{best_l_idx('irpca_lin')                     }.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca')}.iter, results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca')}.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results. rosl.    info{best_l_idx('rosl'     ), best_r_idx('rosl' )}.iter, results. rosl.    info{best_l_idx('rosl'     ), best_r_idx('rosl' )}.err, linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Convergence criterion', 'FontSize', axis_fontsize);
xlabel('Number of iterations', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([1 iter_show]);
ylim([1.0e-7, 1]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_convergence.eps'), 'psc2');

%% ROC curves
cnt = 1;

fig = figure('Position', position); 
fpr = norminv(best_fpr('rpca_alm' )); fnr = norminv(1 - best_tpr('rpca_alm' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
fpr = norminv(best_fpr('rpca_apg' )); fnr = norminv(1 - best_tpr('rpca_apg' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('brpca'    )); fnr = norminv(1 - best_tpr('brpca'    )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('irpca_sub')); fnr = norminv(1 - best_tpr('irpca_sub')); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('irpca_lin')); fnr = norminv(1 - best_tpr('irpca_lin')); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('orpca'    )); fnr = norminv(1 - best_tpr('orpca'    )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('rosl'     )); fnr = norminv(1 - best_tpr('rosl'     )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('False Negative Rate', 'FontSize', axis_fontsize);
xlabel('False Positive Rate', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xtick = [10.^(-4:-1) 0.5]; set(gca, 'XTick', norminv(xtick)); set(gca, 'XTickLabel', {xtick});
ytick = 0.1 : 0.2 : 0.9;   set(gca, 'YTick', norminv(ytick)); set(gca, 'YTickLabel', {ytick});
xlim(norminv([0.0000007 0.76])); 
ylim(norminv([0.05 0.992]));
set(fig, 'PaperPositionMode', 'auto'); 
saveas(fig, fullfile(savepath, 'mat_highway_roc_curves.eps'), 'psc2');

%% best results
fprintf('RPCA (alm): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca_alm')));
fprintf('\t time = %g sec \n', results.rpca.alm.time(best_l_idx('rpca_alm')));
fprintf('\t auc = %g \n', best_auc('rpca_alm'));
fprintf('\t iterations = %d \n\n', results.rpca.alm.info{best_l_idx('rpca_alm')}.iter(end));

fprintf('RPCA (apg): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca_apg')));
fprintf('\t time = %g sec \n', results.rpca.apg.time(best_l_idx('rpca_apg')));
fprintf('\t auc = %g \n', best_auc('rpca_apg'));
fprintf('\t iterations = %d \n\n', results.rpca.apg.info{best_l_idx('rpca_apg')}.iter(end));

fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('brpca')));
fprintf('\t rank = %g \n', mat_rank(best_r_idx('brpca')));
fprintf('\t time = %g sec \n', results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')));
fprintf('\t auc = %g \n', best_auc('brpca'));
fprintf('\t iterations = %d \n\n', results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end));

fprintf('IRPCA (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_sub')));
fprintf('\t time = %g sec \n', results.irpca.sub.time(best_l_idx('irpca_sub')));
fprintf('\t auc = %g \n', best_auc('irpca_sub'));
fprintf('\t iterations = %d \n\n', results.irpca.sub.info{best_l_idx('irpca_sub')}.iter(end));

fprintf('IRPCA (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_lin')));
fprintf('\t time = %g sec \n', results.irpca.lin.time(best_l_idx('irpca_lin')));
fprintf('\t auc = %g \n', best_auc('irpca_lin'));
fprintf('\t iterations = %d \n\n', results.irpca.lin.info{best_l_idx('irpca_lin')}.iter(end));

fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('orpca')));
fprintf('\t rank = %g \n', mat_rank(best_r_idx('orpca')));
fprintf('\t time = %g sec \n', results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')));
fprintf('\t auc = %g \n', best_auc('orpca'));
fprintf('\t iterations = %d \n\n', results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end));

fprintf('ROSL: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rosl')));
fprintf('\t rank = %g \n', mat_rank(best_r_idx('rosl')));
fprintf('\t time = %g sec \n', results.rosl.time(best_l_idx('rosl'),best_r_idx('rosl')));
fprintf('\t auc = %g \n', best_auc('rosl'));
fprintf('\t iterations = %d \n\n', results.rosl.info{best_l_idx('rosl'),best_r_idx('rosl')}.iter(end));

%% print as latex table
fprintf('RPCA (alm) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_auc('rpca_alm'), ...
    results.rpca.alm.time(best_l_idx('rpca_alm')), ...
    results.rpca.alm.info{best_l_idx('rpca_alm')}.iter(end), ...
    lambda(best_l_idx('rpca_alm')) );

fprintf('RPCA (apg) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_auc('rpca_apg'), ...
    results.rpca.apg.time(best_l_idx('rpca_apg')), ...
    results.rpca.apg.info{best_l_idx('rpca_apg')}.iter(end), ...
    lambda(best_l_idx('rpca_apg')) );

fprintf('BRPCA & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_auc('brpca'), ...
    results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')), ...
    results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end), ...
    lambda(best_l_idx('brpca')), ...
    mat_rank(best_r_idx('brpca')) );

fprintf('IRPCA (sub) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_auc('irpca_sub'), ...
    results.irpca.sub.time(best_l_idx('irpca_sub')), ...
    results.irpca.sub.info{best_l_idx('irpca_sub')}.iter(end), ...
    lambda(best_l_idx('irpca_sub')) );

fprintf('IRPCA (lin) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_auc('irpca_lin'), ...
    results.irpca.lin.time(best_l_idx('irpca_lin')), ...
    results.irpca.lin.info{best_l_idx('irpca_lin')}.iter(end), ...
    lambda(best_l_idx('irpca_lin')) );

fprintf('ORPCA & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_auc('orpca'), ...
    results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')), ...
    results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end), ...
    lambda(best_l_idx('orpca')), ...
    mat_rank(best_r_idx('orpca')) );

fprintf('ROSL & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_auc('rosl'), ...
    results.rosl.time(best_l_idx('rosl'),best_r_idx('rosl')), ...
    results.rosl.info{best_l_idx('rosl'),best_r_idx('rosl')}.iter(end), ...
    lambda(best_l_idx('rosl')), ...
    mat_rank(best_r_idx('rosl')) );

%% display backgrounds
disp_imdata(cat(3, X, ...
    results. rpca.alm.A{best_l_idx('rpca_alm' )                     }, ...
    results. rpca.apg.A{best_l_idx('rpca_apg' )                     }, ...
    results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca')}, ...
    results.irpca.sub.A{best_l_idx('irpca_sub')                     }, ...
    results.irpca.lin.A{best_l_idx('irpca_lin')                     }, ...
    results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca')}, ...
    results. rosl.    A{best_l_idx('rosl'     ), best_r_idx('rosl' )} ), ...
    imsize, ...
    {'Original', 'RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL'}, ...
    [2 4], range, inf);

%% display foregrounds
disp_imdata(cat(3, true_mask, ...
    best_foreground('rpca_alm' ), ...
    best_foreground('rpca_apg' ), ...
    best_foreground('brpca'    ), ...
    best_foreground('irpca_sub'), ...
    best_foreground('irpca_lin'), ...
    best_foreground('orpca'    ), ...
    best_foreground('rosl'     ) ), ...
    imsize, ...
    {'Original', 'RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL'}, ...
    [2 4], [], inf);

%% save images
im_idx = [130 164 315 355];
imformat = 'png';

for i = 1:length(im_idx)
    imfile = sprintf('mat_highway_original_%d.%s', im_idx(i), imformat); im = X        (:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_highway_mask_%d.%s'    , im_idx(i), imformat); im = true_mask(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    
    imfile = sprintf('mat_bg_rpca_alm_%d.%s' , im_idx(i), imformat); im = results. rpca.alm.A{best_l_idx('rpca_alm' )                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_bg_rpca_apg_%d.%s' , im_idx(i), imformat); im = results. rpca.apg.A{best_l_idx('rpca_apg' )                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_bg_brpca_%d.%s'    , im_idx(i), imformat); im = results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca')}(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_bg_irpca_sub_%d.%s', im_idx(i), imformat); im = results.irpca.sub.A{best_l_idx('irpca_sub')                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_bg_irpca_lin_%d.%s', im_idx(i), imformat); im = results.irpca.lin.A{best_l_idx('irpca_lin')                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_bg_orpca_%d.%s'    , im_idx(i), imformat); im = results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca')}(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_bg_rosl_%d.%s'     , im_idx(i), imformat); im = results. rosl.    A{best_l_idx('rosl'     ), best_r_idx('rosl' )}(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    
    imfile = sprintf('mat_fg_rpca_alm_%d.%s' , im_idx(i), imformat); im = best_foreground('rpca_alm' ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_fg_rpca_apg_%d.%s' , im_idx(i), imformat); im = best_foreground('rpca_apg' ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_fg_brpca_%d.%s'    , im_idx(i), imformat); im = best_foreground('brpca'    ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_fg_irpca_sub_%d.%s', im_idx(i), imformat); im = best_foreground('irpca_sub'); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_fg_irpca_lin_%d.%s', im_idx(i), imformat); im = best_foreground('irpca_lin'); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_fg_orpca_%d.%s'    , im_idx(i), imformat); im = best_foreground('orpca'    ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_fg_rosl_%d.%s'     , im_idx(i), imformat); im = best_foreground('rosl'     ); im = im(:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
end

%% save movies
movsize = [1 2] .* imsize;
movformat = 'avi';
movfile = sprintf('mat_highway.%s', movformat); mat2movie([X; true_mask], movsize, fullfile(savepath, movfile));

movfile = sprintf('mat_rpca_alm.%s' , movformat); im = best_foreground('rpca_alm' ); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results. rpca.alm.A{best_l_idx('rpca_alm' )                     }; im], movsize, fullfile(savepath, movfile));
movfile = sprintf('mat_rpca_apg.%s' , movformat); im = best_foreground('rpca_apg' ); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results. rpca.apg.A{best_l_idx('rpca_apg' )                     }; im], movsize, fullfile(savepath, movfile));
movfile = sprintf('mat_brpca.%s'    , movformat); im = best_foreground('brpca'    ); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca')}; im], movsize, fullfile(savepath, movfile));
movfile = sprintf('mat_irpca_sub.%s', movformat); im = best_foreground('irpca_sub'); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results.irpca.sub.A{best_l_idx('irpca_sub')                     }; im], movsize, fullfile(savepath, movfile));
movfile = sprintf('mat_irpca_lin.%s', movformat); im = best_foreground('irpca_lin'); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results.irpca.lin.A{best_l_idx('irpca_lin')                     }; im], movsize, fullfile(savepath, movfile));
movfile = sprintf('mat_orpca.%s'    , movformat); im = best_foreground('orpca'    ); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca')}; im], movsize, fullfile(savepath, movfile));
movfile = sprintf('mat_rosl.%s'     , movformat); im = best_foreground('rosl'     ); im = im - min(im(:)); im = im / max(im(:)); mat2movie([results. rosl.    A{best_l_idx('rosl'     ), best_r_idx('rosl' )}; im], movsize, fullfile(savepath, movfile));
