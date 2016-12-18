% plot_results_mat_yale.m
% Plots the results of the denoising experiment on Yale B using matrix
% methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
close all;
clc;

% choose data
dims = [48 42 64];
noise_type = 'no_noise';
noise_param = 40;
method_type = 'matrices';

% load data
path = '/vol/bitbucket/gp1813/experiments/yale';
if strcmp(noise_type, 'no_noise')
    file = sprintf('%dx%dx%d_no_noise.mat', dims);
else
    file = sprintf('%dx%dx%d_%s_%g.mat', dims, noise_type, noise_param);
end
load(fullfile(path, method_type, file));
load(fullfile(path, sprintf('params_%s.mat', noise_type)));

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
load(fullfile(path, 'data', file));
imsize = dims([1 2]);
range = [0 1];

% figure save path
savepath = fullfile('..', 'results', 'yale');

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

% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_err   = containers.Map;
[best_err('rpca_alm' ), best_l_idx('rpca_alm' )] = min(results. rpca.alm.err); 
[best_err('rpca_apg' ), best_l_idx('rpca_apg' )] = min(results. rpca.apg.err); 
[best_err('brpca'    ), best_l_idx('brpca'    )] = min(results.brpca.    err); 
[best_err('irpca_sub'), best_l_idx('irpca_sub')] = min(results.irpca.sub.err); 
[best_err('irpca_lin'), best_l_idx('irpca_lin')] = min(results.irpca.lin.err); 
[best_err('orpca'    ), best_l_idx('orpca'    )] = min(results.orpca.    err); 
[best_err('rosl'     ), best_l_idx('rosl'     )] = min(results. rosl.    err); 

[best_err('brpca'), best_r_idx('brpca')] = min(best_err('brpca')); 
[best_err('orpca'), best_r_idx('orpca')] = min(best_err('orpca')); 
[best_err('rosl' ), best_r_idx('rosl' )] = min(best_err('rosl' )); 

temp = best_l_idx('brpca'); best_l_idx('brpca') = temp(best_r_idx('brpca'));
temp = best_l_idx('orpca'); best_l_idx('orpca') = temp(best_r_idx('orpca'));
temp = best_l_idx('rosl' ); best_l_idx('rosl' ) = temp(best_r_idx('rosl' ));

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
semilogy(results. rpca.alm.info{best_l_idx('rpca_alm' )                     }.iter, results. rpca.alm.info{best_l_idx('rpca_alm' )                     }.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
semilogy(results. rpca.apg.info{best_l_idx('rpca_apg' )                     }.iter, results. rpca.apg.info{best_l_idx('rpca_apg' )                     }.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca')}.iter, results.brpca.    info{best_l_idx('brpca'    ), best_r_idx('brpca')}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.sub.info{best_l_idx('irpca_sub')                     }.iter, results.irpca.sub.info{best_l_idx('irpca_sub')                     }.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.irpca.lin.info{best_l_idx('irpca_lin')                     }.iter, results.irpca.lin.info{best_l_idx('irpca_lin')                     }.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca')}.iter, results.orpca.    info{best_l_idx('orpca'    ), best_r_idx('orpca')}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
semilogy(results. rosl.    info{best_l_idx('rosl'     ), best_r_idx('rosl' )}.iter, results. rosl.    info{best_l_idx('rosl'     ), best_r_idx('rosl' )}.err,  linestyle, 'Marker', marker{cnt}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
ylabel('Convergence criterion', 'FontSize', axis_fontsize);
xlabel('Number of iterations', 'FontSize', axis_fontsize);
legend('RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL', 'Location', 'NorthEastOutside');
set(gca, 'FontSize', axis_fontsize);
xlim([1 iter_show]);
ylim([1.0e-7, 1]);
set(fig, 'PaperPositionMode', 'auto'); 
%saveas(fig, fullfile(savepath, 'synthetic_mat_convergence.eps'), 'psc2');

%% best results
fprintf('RPCA (alm): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca_alm')));
fprintf('\t time = %g sec \n', results.rpca.alm.time(best_l_idx('rpca_alm')));
fprintf('\t error = %g \n', best_err('rpca_alm'));
fprintf('\t iterations = %d \n\n', results.rpca.alm.info{best_l_idx('rpca_alm')}.iter(end));

fprintf('RPCA (apg): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca_apg')));
fprintf('\t time = %g sec \n', results.rpca.apg.time(best_l_idx('rpca_apg')));
fprintf('\t error = %g \n', best_err('rpca_apg'));
fprintf('\t iterations = %d \n\n', results.rpca.apg.info{best_l_idx('rpca_apg')}.iter(end));

fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('brpca')));
fprintf('\t rank = %g \n', mat_rank(best_r_idx('brpca')));
fprintf('\t time = %g sec \n', results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')));
fprintf('\t error = %g \n', best_err('brpca'));
fprintf('\t iterations = %d \n\n', results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end));

fprintf('IRPCA (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_sub')));
fprintf('\t time = %g sec \n', results.irpca.sub.time(best_l_idx('irpca_sub')));
fprintf('\t error = %g \n', best_err('irpca_sub'));
fprintf('\t iterations = %d \n\n', results.irpca.sub.info{best_l_idx('irpca_sub')}.iter(end));

fprintf('IRPCA (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_lin')));
fprintf('\t time = %g sec \n', results.irpca.lin.time(best_l_idx('irpca_lin')));
fprintf('\t error = %g \n', best_err('irpca_lin'));
fprintf('\t iterations = %d \n\n', results.irpca.lin.info{best_l_idx('irpca_lin')}.iter(end));

fprintf('ORPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('orpca')));
fprintf('\t rank = %g \n', mat_rank(best_r_idx('orpca')));
fprintf('\t time = %g sec \n', results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')));
fprintf('\t error = %g \n', best_err('orpca'));
fprintf('\t iterations = %d \n\n', results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end));

fprintf('ROSL: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rosl')));
fprintf('\t rank = %g \n', mat_rank(best_r_idx('rosl')));
fprintf('\t time = %g sec \n', results.rosl.time(best_l_idx('rosl'),best_r_idx('rosl')));
fprintf('\t error = %g \n', best_err('rosl'));
fprintf('\t iterations = %d \n\n', results.rosl.info{best_l_idx('rosl'),best_r_idx('rosl')}.iter(end));

%% print as latex table
fprintf('RPCA (alm) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_err('rpca_alm'), ...
    results.rpca.alm.time(best_l_idx('rpca_alm')), ...
    results.rpca.alm.info{best_l_idx('rpca_alm')}.iter(end), ...
    lambda(best_l_idx('rpca_alm')) );

fprintf('RPCA (apg) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_err('rpca_apg'), ...
    results.rpca.apg.time(best_l_idx('rpca_apg')), ...
    results.rpca.apg.info{best_l_idx('rpca_apg')}.iter(end), ...
    lambda(best_l_idx('rpca_apg')) );

fprintf('BRPCA & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_err('brpca'), ...
    results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')), ...
    results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end), ...
    lambda(best_l_idx('brpca')), ...
    mat_rank(best_r_idx('brpca')) );

fprintf('IRPCA (sub) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_err('irpca_sub'), ...
    results.irpca.sub.time(best_l_idx('irpca_sub')), ...
    results.irpca.sub.info{best_l_idx('irpca_sub')}.iter(end), ...
    lambda(best_l_idx('irpca_sub')) );

fprintf('IRPCA (lin) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_err('irpca_lin'), ...
    results.irpca.lin.time(best_l_idx('irpca_lin')), ...
    results.irpca.lin.info{best_l_idx('irpca_lin')}.iter(end), ...
    lambda(best_l_idx('irpca_lin')) );

fprintf('ORPCA & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_err('orpca'), ...
    results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')), ...
    results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end), ...
    lambda(best_l_idx('orpca')), ...
    mat_rank(best_r_idx('orpca')) );

fprintf('ROSL & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_err('rosl'), ...
    results.rosl.time(best_l_idx('rosl'),best_r_idx('rosl')), ...
    results.rosl.info{best_l_idx('rosl'),best_r_idx('rosl')}.iter(end), ...
    lambda(best_l_idx('rosl')), ...
    mat_rank(best_r_idx('rosl')) );

%% display images
disp_imdata(cat(3, arr2mat(thisX_noise, [1 2]), ...
    results. rpca.alm.A{best_l_idx('rpca_alm' )                     }, ...
    results. rpca.apg.A{best_l_idx('rpca_apg' )                     }, ...
    results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca')}, ...
    results.irpca.sub.A{best_l_idx('irpca_sub')                     }, ...
    results.irpca.lin.A{best_l_idx('irpca_lin')                     }, ...
    results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca')}, ...
    results. rosl.    A{best_l_idx('rosl'     ), best_r_idx('rosl' )} ), ...
    imsize, ...
    {'Original', 'RPCA (alm)', 'RPCA (apg)', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'ORPCA', 'ROSL'}, ...
    [2 4], range);

%% save images
im_idx = [1 7 30];
imformat = 'png';

if length(im_idx) == 1
    imfile = sprintf('mat_%s_%g_original.%s' , noise_type, noise_param, imformat); im = thisX_noise                                                    (:,:,im_idx); imwrite(        im,          fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_rpca_alm.%s' , noise_type, noise_param, imformat); im = results. rpca.alm.A{best_l_idx('rpca_alm' )                     }(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_rpca_apg.%s' , noise_type, noise_param, imformat); im = results. rpca.apg.A{best_l_idx('rpca_apg' )                     }(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_brpca.%s'    , noise_type, noise_param, imformat); im = results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca')}(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_irpca_sub.%s', noise_type, noise_param, imformat); im = results.irpca.sub.A{best_l_idx('irpca_sub')                     }(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_irpca_lin.%s', noise_type, noise_param, imformat); im = results.irpca.lin.A{best_l_idx('irpca_lin')                     }(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_orpca.%s'    , noise_type, noise_param, imformat); im = results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca')}(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    imfile = sprintf('mat_%s_%g_rosl.%s'     , noise_type, noise_param, imformat); im = results. rosl.    A{best_l_idx('rosl'     ), best_r_idx('rosl' )}(:,im_idx); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
else
    for i = 1:length(im_idx)
        imfile = sprintf('mat_%s_%g_original_%d.%s' , noise_type, noise_param, im_idx(i), imformat); im = thisX_noise                                                    (:,:,im_idx(i)); imwrite(        im,          fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_rpca_alm_%d.%s' , noise_type, noise_param, im_idx(i), imformat); im = results. rpca.alm.A{best_l_idx('rpca_alm' )                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_rpca_apg_%d.%s' , noise_type, noise_param, im_idx(i), imformat); im = results. rpca.apg.A{best_l_idx('rpca_apg' )                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_brpca_%d.%s'    , noise_type, noise_param, im_idx(i), imformat); im = results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca')}(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_irpca_sub_%d.%s', noise_type, noise_param, im_idx(i), imformat); im = results.irpca.sub.A{best_l_idx('irpca_sub')                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_irpca_lin_%d.%s', noise_type, noise_param, im_idx(i), imformat); im = results.irpca.lin.A{best_l_idx('irpca_lin')                     }(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_orpca_%d.%s'    , noise_type, noise_param, im_idx(i), imformat); im = results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca')}(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
        imfile = sprintf('mat_%s_%g_rosl_%d.%s'     , noise_type, noise_param, im_idx(i), imformat); im = results. rosl.    A{best_l_idx('rosl'     ), best_r_idx('rosl' )}(:,im_idx(i)); imwrite(reshape(im, imsize), fullfile(savepath, imfile), imformat);
    end
end
