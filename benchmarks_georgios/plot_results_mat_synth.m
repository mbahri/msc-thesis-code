% plot_results_mat_synth.m
% Plots the results of the experiment on synthetic data using matrix
% methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
close all;
clc;

% choose data
dims = [1000 1000];
rank = 0.05;
data_magn = 1;
err_perc = 0.10;
err_magn = 100;

% load data
path = '/vol/bitbucket/gp1813/experiments/synthetic/matrices';
file = sprintf('%dx%d_rank_%g_data_%g_err_%g_%d_v2.mat', dims, rank, data_magn, err_perc, err_magn);
load(fullfile(path, file));
load(fullfile(path, 'params_v2.mat'));

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

% figure save path
savepath = fullfile('..', 'results', 'synthetic');

% p parameter
p_idx = 3;

% best results for each method
best_l_idx = containers.Map;
best_err   = containers.Map;
[best_err('rpca_alm' ), best_l_idx('rpca_alm' )] = min(results. rpca.alm.err(:,p_idx)); 
[best_err('rpca_apg' ), best_l_idx('rpca_apg' )] = min(results. rpca.apg.err(:,p_idx)); 
[best_err('brpca'    ), best_l_idx('brpca'    )] = min(results.brpca.    err(:,p_idx)); 
[best_err('irpca_sub'), best_l_idx('irpca_sub')] = min(results.irpca.sub.err(:,p_idx)); 
[best_err('irpca_lin'), best_l_idx('irpca_lin')] = min(results.irpca.lin.err(:,p_idx)); 
[best_err('orpca'    ), best_l_idx('orpca'    )] = min(results.orpca.    err(:,p_idx)); 
[best_err('rosl'     ), best_l_idx('rosl'     )] = min(results. rosl.    err(:,p_idx)); 
maxiter = max([results. rpca.alm.info{best_l_idx('rpca_alm' ),p_idx}.iter(end) ...
               results. rpca.apg.info{best_l_idx('rpca_apg' ),p_idx}.iter(end) ...
               results.brpca.    info{best_l_idx('brpca'    ),p_idx}.iter(end) ...
               results.irpca.sub.info{best_l_idx('irpca_sub'),p_idx}.iter(end) ...
               results.irpca.lin.info{best_l_idx('irpca_lin'),p_idx}.iter(end) ...
               results.orpca.    info{best_l_idx('orpca'    ),p_idx}.iter(end) ...
               results. rosl.    info{best_l_idx('rosl'     ),p_idx}.iter(end)] );

%% time plots
cnt = 1;

fig = figure('Position', position);
semilogx(lambda, results. rpca.alm.time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogx(lambda, results. rpca.apg.time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.brpca.    time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.sub.time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.irpca.lin.time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results.orpca.    time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogx(lambda, results. rosl.    time(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Total execution time [sec]', 'FontSize', axis_fontsize);
xlabel('\lambda', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([min(lambda) max(lambda)]);
set(fig, 'PaperPositionMode', 'auto'); 
saveas(fig, fullfile(savepath, 'synthetic_mat_time.eps'), 'psc2');

%% performance plots
cnt = 1;

fig = figure('Position', position);
loglog(lambda, results. rpca.alm.err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
loglog(lambda, results. rpca.apg.err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.brpca.    err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.irpca.sub.err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.irpca.lin.err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results.orpca.    err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
loglog(lambda, results. rosl.    err(:,p_idx), linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_big, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Relative recovery error',   'FontSize', axis_fontsize);
xlabel('\lambda', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([min(lambda) max(lambda)]);
ylim([1.0e-7 1.0e+2]);
set(fig, 'PaperPositionMode', 'auto'); 
saveas(fig, fullfile(savepath, 'synthetic_mat_performance.eps'), 'psc2');

%% convergence plots
iter_show = maxiter; cnt = 1;

fig = figure('Position', position); 
semilogy(results. rpca.alm.info{best_l_idx('rpca_alm' ),p_idx}.iter, results. rpca.alm.info{best_l_idx('rpca_alm' ),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogy(results. rpca.apg.info{best_l_idx('rpca_apg' ),p_idx}.iter, results. rpca.apg.info{best_l_idx('rpca_apg' ),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.brpca.    info{best_l_idx('brpca'    ),p_idx}.iter, results.brpca.    info{best_l_idx('brpca'    ),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.sub.info{best_l_idx('irpca_sub'),p_idx}.iter, results.irpca.sub.info{best_l_idx('irpca_sub'),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.lin.info{best_l_idx('irpca_lin'),p_idx}.iter, results.irpca.lin.info{best_l_idx('irpca_lin'),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.orpca.    info{best_l_idx('orpca'    ),p_idx}.iter, results.orpca.    info{best_l_idx('orpca'    ),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results. rosl.    info{best_l_idx('rosl'     ),p_idx}.iter, results. rosl.    info{best_l_idx('rosl'     ),p_idx}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Convergence criterion', 'FontSize', axis_fontsize);
xlabel('Number of iterations', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([1 iter_show]);
ylim([1.0e-7, 1]);
set(fig, 'PaperPositionMode', 'auto'); 
saveas(fig, fullfile(savepath, 'synthetic_mat_convergence.eps'), 'psc2');

%% best results
fprintf('RPCA (alm): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca_alm')));
fprintf('\t time = %g sec \n', results.rpca.alm.time(best_l_idx('rpca_alm'),p_idx));
fprintf('\t error = %g \n', best_err('rpca_alm'));
fprintf('\t iterations = %d \n\n', results.rpca.alm.info{best_l_idx('rpca_alm'),p_idx}.iter(end));

fprintf('RPCA (apg): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca_apg')));
fprintf('\t time = %g sec \n', results.rpca.apg.time(best_l_idx('rpca_apg'),p_idx));
fprintf('\t error = %g \n', best_err('rpca_apg'));
fprintf('\t iterations = %d \n\n', results.rpca.apg.info{best_l_idx('rpca_apg'),p_idx}.iter(end));

fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('brpca')));
fprintf('\t time = %g sec \n', results.brpca.time(best_l_idx('brpca'),p_idx));
fprintf('\t error = %g \n', best_err('brpca'));
fprintf('\t iterations = %d \n\n', results.brpca.info{best_l_idx('brpca'),p_idx}.iter(end));

fprintf('IRPCA (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_sub')));
fprintf('\t time = %g sec \n', results.irpca.sub.time(best_l_idx('irpca_sub'),p_idx));
fprintf('\t error = %g \n', best_err('irpca_sub'));
fprintf('\t iterations = %d \n\n', results.irpca.sub.info{best_l_idx('irpca_sub'),p_idx}.iter(end));

fprintf('IRPCA (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_lin')));
fprintf('\t time = %g sec \n', results.irpca.lin.time(best_l_idx('irpca_lin'),p_idx));
fprintf('\t error = %g \n', best_err('irpca_lin'));
fprintf('\t iterations = %d \n\n', results.irpca.lin.info{best_l_idx('irpca_lin'),p_idx}.iter(end));

fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('orpca')));
fprintf('\t time = %g sec \n', results.orpca.time(best_l_idx('orpca'),p_idx));
fprintf('\t error = %g \n', best_err('orpca'));
fprintf('\t iterations = %d \n\n', results.orpca.info{best_l_idx('orpca'),p_idx}.iter(end));

fprintf('ROSL: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rosl')));
fprintf('\t time = %g sec \n', results.rosl.time(best_l_idx('rosl'),p_idx));
fprintf('\t error = %g \n', best_err('rosl'));
fprintf('\t iterations = %d \n\n', results.rosl.info{best_l_idx('rosl'),p_idx}.iter(end));

%% print as latex table
fprintf('RPCA (alm) & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('rpca_alm'), ...
    results.rpca.alm.time(best_l_idx('rpca_alm'),p_idx), ...
    results.rpca.alm.info{best_l_idx('rpca_alm'),p_idx}.iter(end), ...
    lambda(best_l_idx('rpca_alm')) );

fprintf('RPCA (apg) & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('rpca_apg'), ...
    results.rpca.apg.time(best_l_idx('rpca_apg'),p_idx), ...
    results.rpca.apg.info{best_l_idx('rpca_apg'),p_idx}.iter(end), ...
    lambda(best_l_idx('rpca_apg')) );

fprintf('BRPCA & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('brpca'), ...
    results.brpca.time(best_l_idx('brpca'),p_idx), ...
    results.brpca.info{best_l_idx('brpca'),p_idx}.iter(end), ...
    lambda(best_l_idx('brpca')) );

fprintf('IRPCA (sub) & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('irpca_sub'), ...
    results.irpca.sub.time(best_l_idx('irpca_sub'),p_idx), ...
    results.irpca.sub.info{best_l_idx('irpca_sub'),p_idx}.iter(end), ...
    lambda(best_l_idx('irpca_sub')) );

fprintf('IRPCA (lin) & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('irpca_lin'), ...
    results.irpca.lin.time(best_l_idx('irpca_lin'),p_idx), ...
    results.irpca.lin.info{best_l_idx('irpca_lin'),p_idx}.iter(end), ...
    lambda(best_l_idx('irpca_lin')) );

fprintf('ORPCA & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('orpca'), ...
    results.orpca.time(best_l_idx('orpca'),p_idx), ...
    results.orpca.info{best_l_idx('orpca'),p_idx}.iter(end), ...
    lambda(best_l_idx('orpca')) );

fprintf('ROSL & $%g$ & $%.3f$ & $%d$ & $%.4f$ \\\\ \n', ...
    best_err('rosl'), ...
    results.rosl.time(best_l_idx('rosl'),p_idx), ...
    results.rosl.info{best_l_idx('rosl'),p_idx}.iter(end), ...
    lambda(best_l_idx('rosl')) );
