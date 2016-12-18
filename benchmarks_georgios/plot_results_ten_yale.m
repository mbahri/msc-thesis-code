% plot_results_ten_yale.m
% Plots the results of the denoising experiment on Yale B using tensor 
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
noise_type = 'random_patch';
noise_param = 40;
method_type = 'tensors';

% load data
path = '/vol/bitbucket/gp1813/experiments/yale';
if length(dims) == 3
    if strcmp(noise_type, 'no_noise')
        file = sprintf('%dx%dx%d_no_noise.mat', dims);
    else
        file = sprintf('%dx%dx%d_%s_%g.mat', dims, noise_type, noise_param);
    end
else
    if strcmp(noise_type, 'no_noise')
        file = sprintf('%dx%dx%dx%d_no_noise.mat', dims);
    else
        file = sprintf('%dx%dx%dx%d_%s_%g.mat', dims, noise_type, noise_param);
    end
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

% display pics parameters
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

% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_err   = containers.Map;
[best_err('rpca'     ), best_l_idx('rpca'     )] = min(results. rpca.    err); 
[best_err('brpca'    ), best_l_idx('brpca'    )] = min(results.brpca.    err); 
[best_err('irpca_sub'), best_l_idx('irpca_sub')] = min(results.irpca.sub.err); 
[best_err('irpca_lin'), best_l_idx('irpca_lin')] = min(results.irpca.lin.err); 
[best_err('orpca'    ), best_l_idx('orpca'    )] = min(results.orpca.    err); 
[best_err('rcpd_sub' ), best_l_idx('rcpd_sub' )] = min(results. rcpd.sub.err); 
[best_err('rcpd_lin' ), best_l_idx('rcpd_lin' )] = min(results. rcpd.lin.err); 

[best_err('brpca'    ), best_r_idx('brpca'    )] = min(best_err('brpca'    )); 
[best_err('irpca_sub'), best_r_idx('irpca_sub')] = min(best_err('irpca_sub')); 
[best_err('irpca_lin'), best_r_idx('irpca_lin')] = min(best_err('irpca_lin')); 
[best_err('orpca'    ), best_r_idx('orpca'    )] = min(best_err('orpca'    )); 
[best_err('rcpd_sub' ), best_r_idx('rcpd_sub' )] = min(best_err('rcpd_sub' )); 
[best_err('rcpd_lin' ), best_r_idx('rcpd_lin' )] = min(best_err('rcpd_lin' )); 

temp = best_l_idx('brpca'    ); best_l_idx('brpca'    ) = temp(best_r_idx('brpca'    ));
temp = best_l_idx('irpca_sub'); best_l_idx('irpca_sub') = temp(best_r_idx('irpca_sub'));
temp = best_l_idx('irpca_lin'); best_l_idx('irpca_lin') = temp(best_r_idx('irpca_lin'));
temp = best_l_idx('orpca'    ); best_l_idx('orpca'    ) = temp(best_r_idx('orpca'    ));
temp = best_l_idx('rcpd_sub' ); best_l_idx('rcpd_sub' ) = temp(best_r_idx('rcpd_sub' ));
temp = best_l_idx('rcpd_lin' ); best_l_idx('rcpd_lin' ) = temp(best_r_idx('rcpd_lin' ));

maxiter = max([results. rpca.    info{best_l_idx('rpca'     )                        }.iter(end) ...
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
fprintf('\t error = %g \n', best_err('rpca'));
fprintf('\t iterations = %d \n\n', results.rpca.info{best_l_idx('rpca')}.iter(end));

fprintf('BRPCA: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('brpca')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('brpca')));
fprintf('\t time = %g sec \n', results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')));
fprintf('\t error = %g \n', best_err('brpca'));
fprintf('\t iterations = %d \n\n', results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end));

fprintf('IRPCA (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_sub')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('irpca_sub')));
fprintf('\t time = %g sec \n', results.irpca.sub.time(best_l_idx('irpca_sub'),best_r_idx('irpca_sub')));
fprintf('\t error = %g \n', best_err('irpca_sub'));
fprintf('\t iterations = %d \n\n', results.irpca.sub.info{best_l_idx('irpca_sub'),best_r_idx('irpca_sub')}.iter(end));

fprintf('IRPCA (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_lin')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('irpca_lin')));
fprintf('\t time = %g sec \n', results.irpca.lin.time(best_l_idx('irpca_lin'),best_r_idx('irpca_lin')));
fprintf('\t error = %g \n', best_err('irpca_lin'));
fprintf('\t iterations = %d \n\n', results.irpca.lin.info{best_l_idx('irpca_lin'),best_r_idx('irpca_lin')}.iter(end));

fprintf('RHOSVD: \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('orpca')));
fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('orpca')));
fprintf('\t time = %g sec \n', results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')));
fprintf('\t error = %g \n', best_err('orpca'));
fprintf('\t iterations = %d \n\n', results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end));

fprintf('RCPD (sub): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rcpd_sub')));
fprintf('\t rank = %g \n', ten_rank(best_r_idx('rcpd_sub')));
fprintf('\t time = %g sec \n', results.rcpd.sub.time(best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')));
fprintf('\t error = %g \n', best_err('rcpd_sub'));
fprintf('\t iterations = %d \n\n', results.rcpd.sub.info{best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')}.iter(end));

fprintf('RCPD (lin): \n');
fprintf('\t lambda = %g \n', lambda(best_l_idx('rcpd_lin')));
fprintf('\t rank = %g \n', ten_rank(best_r_idx('rcpd_lin')));
fprintf('\t time = %g sec \n', results.rcpd.lin.time(best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')));
fprintf('\t error = %g \n', best_err('rcpd_lin'));
fprintf('\t iterations = %d \n\n', results.rcpd.lin.info{best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')}.iter(end));

%% print as latex table
fprintf('RPCA & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & -- \\\\ \n', ...
    best_err('rpca'), ...
    results.rpca.time(best_l_idx('rpca')), ...
    results.rpca.info{best_l_idx('rpca')}.iter(end), ...
    lambda(best_l_idx('rpca')) );

fprintf('BRPCA & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%g\\times$ \\\\ \n', ...
    best_err('brpca'), ...
    results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')), ...
    results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end), ...
    lambda(best_l_idx('brpca')), ...
    ten_nrank_perc(best_r_idx('brpca')) );

fprintf('IRPCA (sub) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%g\\times$ \\\\ \n', ...
    best_err('irpca_sub'), ...
    results.irpca.sub.time(best_l_idx('irpca_sub'),best_r_idx('irpca_sub')), ...
    results.irpca.sub.info{best_l_idx('irpca_sub'),best_r_idx('irpca_sub')}.iter(end), ...
    lambda(best_l_idx('irpca_sub')), ...
    ten_nrank_perc(best_r_idx('irpca_sub')) );

fprintf('IRPCA (lin) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%g\\times$ \\\\ \n', ...
    best_err('irpca_lin'), ...
    results.irpca.lin.time(best_l_idx('irpca_lin'),best_r_idx('irpca_lin')), ...
    results.irpca.lin.info{best_l_idx('irpca_lin'),best_r_idx('irpca_lin')}.iter(end), ...
    lambda(best_l_idx('irpca_lin')), ...
    ten_nrank_perc(best_r_idx('irpca_lin')) );

fprintf('RHOSVD & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%g\\times$ \\\\ \n', ...
    best_err('orpca'), ...
    results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')), ...
    results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end), ...
    lambda(best_l_idx('orpca')), ...
    ten_nrank_perc(best_r_idx('orpca')) );

fprintf('RCPD (sub) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_err('rcpd_sub'), ...
    results.rcpd.sub.time(best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')), ...
    results.rcpd.sub.info{best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')}.iter(end), ...
    lambda(best_l_idx('rcpd_sub')), ...
    ten_rank(best_r_idx('rcpd_sub')) );

fprintf('RCPD (lin) & $%.3f$ & $%.3f$ & $%d$ & $%.4f$ & $%d$ \\\\ \n', ...
    best_err('rcpd_lin'), ...
    results.rcpd.lin.time(best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')), ...
    results.rcpd.lin.info{best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')}.iter(end), ...
    lambda(best_l_idx('rcpd_lin')), ...
    ten_rank(best_r_idx('rcpd_lin')) );

%% display images
disp_imdata(cat(3, arr2mat(thisX_noise, [1 2]), ...
    arr2mat(results. rpca.    A{best_l_idx('rpca'     )                         }, [1 2]), ...
    arr2mat(results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca'    )}, [1 2]), ...
    arr2mat(results.irpca.sub.A{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}, [1 2]), ...
    arr2mat(results.irpca.lin.A{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}, [1 2]), ...
    arr2mat(results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca'    )}, [1 2]), ...
    arr2mat(results. rcpd.sub.A{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}, [1 2]), ...
    arr2mat(results. rcpd.lin.A{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}, [1 2]) ), ...
    imsize, ...
    {'Original', 'RPCA', 'BRPCA', 'IRPCA (sub)', 'IRPCA (lin)', 'RHOSVD', 'RCPD (sub)', 'RCPD (lin)'}, ...
    [2 4], range);

%% save images
im_idx = [1 7 30];
imformat = 'png';

if length(im_idx) == 1
    imfile = sprintf('ten_%s_%g_original.%s' , noise_type, noise_param, imformat); im = thisX_noise                                                          (:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_rpca.%s'     , noise_type, noise_param, imformat); im = results. rpca.    A{best_l_idx('rpca'     )                         }(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_brpca.%s'    , noise_type, noise_param, imformat); im = results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca'    )}(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_irpca_sub.%s', noise_type, noise_param, imformat); im = results.irpca.sub.A{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_irpca_lin.%s', noise_type, noise_param, imformat); im = results.irpca.lin.A{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_orpca.%s'    , noise_type, noise_param, imformat); im = results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca'    )}(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_rcpd_sub.%s' , noise_type, noise_param, imformat); im = results. rcpd.sub.A{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
    imfile = sprintf('ten_%s_%g_rcpd_lin.%s' , noise_type, noise_param, imformat); im = results. rcpd.lin.A{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}(:,:,im_idx); imwrite(double(im), fullfile(savepath, imfile), imformat);
else
    for i = 1:length(im_idx)
        imfile = sprintf('ten_%s_%g_original_%d.%s' , noise_type, noise_param, im_idx(i), imformat); im = thisX_noise                                                          (:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_rpca_%d.%s'     , noise_type, noise_param, im_idx(i), imformat); im = results. rpca.    A{best_l_idx('rpca'     )                         }(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_brpca_%d.%s'    , noise_type, noise_param, im_idx(i), imformat); im = results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca'    )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_irpca_sub_%d.%s', noise_type, noise_param, im_idx(i), imformat); im = results.irpca.sub.A{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_irpca_lin_%d.%s', noise_type, noise_param, im_idx(i), imformat); im = results.irpca.lin.A{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_orpca_%d.%s'    , noise_type, noise_param, im_idx(i), imformat); im = results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca'    )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_rcpd_sub_%d.%s' , noise_type, noise_param, im_idx(i), imformat); im = results. rcpd.sub.A{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
        imfile = sprintf('ten_%s_%g_rcpd_lin_%d.%s' , noise_type, noise_param, im_idx(i), imformat); im = results. rcpd.lin.A{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
    end
end