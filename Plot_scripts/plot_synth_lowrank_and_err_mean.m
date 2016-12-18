% [fg, tpr, fpr, auc, l_idx, r_idx] = select_best_foreground({D.L(:,:,50)}, O(:,:,50), true_mask, linspace(0, 1, 101));
% 
% elastic net: take weighted sum (convex) of l1 and l2 solutions?

%% Initialize and load the parameters
clear
close all

legend_style = @(LineStyle, LineWidth, MarkerStyle) plot(0,0,LineStyle,'LineWidth',LineWidth,'Marker',MarkerStyle,'visible','off');


load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_with_mean_rpca2d_gl_l1_sp_1.000000.mat')
[results_gl_l1, best_gl_l1] = distilate(results_mehdi.rpca2d_gl_l1, 'lambda');
[results_gl_l1_err, best_gl_l1_err] = distilate_err(results_mehdi.rpca2d_gl_l1, 'lambda');
[results_gl_l1_mean, best_gl_l1_mean] = distilate_mean(results_mehdi.rpca2d_gl_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_with_mean_rpca2d_gl_l2_sp_1.000000.mat')
[results_gl_l2, best_gl_l2] = distilate(results_mehdi.rpca2d_gl_l2, 'lambda');
[results_gl_l2_err, best_gl_l2_err] = distilate_err(results_mehdi.rpca2d_gl_l2, 'lambda');
[results_gl_l2_mean, best_gl_l2_mean] = distilate_mean(results_mehdi.rpca2d_gl_l2, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_with_mean_rpca2d_l1_sp_1.000000.mat')
[results_l1, best_l1] = distilate(results_mehdi.rpca2d_l1, 'lambda');
[results_l1_err, best_l1_err] = distilate_err(results_mehdi.rpca2d_l1, 'lambda');
[results_l1_mean, best_l1_mean] = distilate_mean(results_mehdi.rpca2d_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_with_mean_rpca2d_l2_sp_1.000000.mat')
[results_l2, best_l2] = distilate(results_mehdi.rpca2d_l2, 'lambda');
[results_l2_err, best_l2_err] = distilate_err(results_mehdi.rpca2d_l2, 'lambda');
[results_l2_mean, best_l2_mean] = distilate_mean(results_mehdi.rpca2d_l2, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_with_mean_rpca2d_gl_l1_sp_1.000000.mat')
[results_gl_l1_60, best_gl_l1_60] = distilate(results_mehdi.rpca2d_gl_l1, 'lambda');
[results_gl_l1_60_err, best_gl_l1_60_err] = distilate_err(results_mehdi.rpca2d_gl_l1, 'lambda');
[results_gl_l1_60_mean, best_gl_l1_60_mean] = distilate_mean(results_mehdi.rpca2d_gl_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_with_mean_rpca2d_gl_l2_sp_1.000000.mat')
[results_gl_l2_60, best_gl_l2_60] = distilate(results_mehdi.rpca2d_gl_l2, 'lambda');
[results_gl_l2_60_err, best_gl_l2_60_err] = distilate_err(results_mehdi.rpca2d_gl_l2, 'lambda');
[results_gl_l2_60_mean, best_gl_l2_60_mean] = distilate_mean(results_mehdi.rpca2d_gl_l2, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_with_mean_rpca2d_l1_sp_1.000000.mat')
[results_l1_60, best_l1_60] = distilate(results_mehdi.rpca2d_l1, 'lambda');
[results_l1_60_err, best_l1_60_err] = distilate_err(results_mehdi.rpca2d_l1, 'lambda');
[results_l1_60_mean, best_l1_60_mean] = distilate_mean(results_mehdi.rpca2d_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_with_mean_rpca2d_l2_sp_1.000000.mat')
[results_l2_60, best_l2_60] = distilate(results_mehdi.rpca2d_l2, 'lambda');
[results_l2_60_err, best_l2_60_err] = distilate_err(results_mehdi.rpca2d_l2, 'lambda');
[results_l2_60_mean, best_l2_60_mean] = distilate_mean(results_mehdi.rpca2d_l2, 'lambda');


[O30, X30] = synthetic_data_30(1000,1000,100,23,41);
[O60, X60] = synthetic_data_60(1000,1000,100,23,41);


%% Plot PSNR and relative error for 30% noise
fig_err_psnr_30 = figure;

subplot(1,2,1);
gca_psnr_30 = gca;
loglog(results_gl_l1.param1, results_gl_l1.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2.param1, results_gl_l2.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1.param1,    results_l1.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2.param1,    results_l2.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_psnr_30.XTick = [1e-3, 1e-2, 1e-1];
gca_psnr_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 30%
ylim([5e-8, 1])
gca_psnr_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_psnr_30.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('\lambda'); ylabel('Rel. Err. - L (30% noise)')

subplot(1,2,2)
gca_psnr_err_30 = gca;
loglog(results_gl_l1.param1, results_gl_l1_err.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2.param1, results_gl_l2_err.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1.param1,    results_l1_err.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2.param1,    results_l2_err.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_psnr_err_30.XTick = [1e-3, 1e-2, 1e-1];
gca_psnr_err_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 30%
ylim([1e-8, 1])
gca_psnr_err_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_psnr_err_30.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('\lambda'); ylabel('Rel. Err. - E (30% noise)')

%tightfig; 

fig_err_psnr_30.Units = 'centimeters';
fig_err_psnr_30.Position = [0 0 15.9 5.6];
fig_err_psnr_30.PaperUnits = 'centimeters';
fig_err_psnr_30.PaperPosition = [0 0 16 6];

saveas(fig_err_psnr_30, 'synth_30_error_mean', 'epsc');

%% Plot PSNR and relative error for 60% noise
fig_err_psnr_60 = figure;

subplot(1,2,1);
gca_psnr_60 = gca;
loglog(results_gl_l1_60.param1, results_gl_l1_60.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2_60.param1, results_gl_l2_60.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1_60.param1,    results_l1_60.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2_60.param1,    results_l2_60.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_psnr_60.XTick = [1e-3, 1e-2, 1e-1];
gca_psnr_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 60%
ylim([1e-2, 1]);
gca_psnr_60.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_psnr_60.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('\lambda'); ylabel('Rel. Err. - L (60% noise)')

subplot(1,2,2)
gca_err_60 = gca;
loglog(results_gl_l1_60.param1, results_gl_l1_60_err.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2_60.param1, results_gl_l2_60_err.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1_60.param1,    results_l1_60_err.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2_60.param1,    results_l2_60_err.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_err_60.XTick = [1e-3, 1e-2, 1e-1];
gca_err_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 60%
ylim([5e-3, 1]);
gca_err_60.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_err_60.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('\lambda'); ylabel('Rel. Err. - E (60% noise)')

%tightfig; 

fig_err_psnr_60.Units = 'centimeters';
fig_err_psnr_60.Position = [0 0 15.9 5.6];
fig_err_psnr_60.PaperUnits = 'centimeters';
fig_err_psnr_60.PaperPosition = [0 0 16 6];

saveas(fig_err_psnr_60, 'synth_60_error_mean', 'epsc');

%% FSIM and sparsity for 30% noise

fig_structure_30 = figure;

% subplot(1,2,1);
% gca_ssim_30 = gca;
% semilogx(results_gl_l1.param1, results_gl_l1.ssim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
% semilogx(results_gl_l2.param1, results_gl_l2.ssim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
% semilogx(results_l1.param1,    results_l1.ssim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
% semilogx(results_l2.param1,    results_l2.ssim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);
% 
% legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')
% 
% xlim([2e-4, 2e-1])
% gca_ssim_30.XTick = [1e-3, 1e-2, 1e-1];
% gca_ssim_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};
% 
% xlabel('\lambda'); ylabel('SSIM (30% noise)')

subplot(1,2,1);
gca_fsim_30 = gca;
semilogx(results_gl_l1.param1, results_gl_l1.fsim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2.param1, results_gl_l2.fsim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1.param1,    results_l1.fsim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2.param1,    results_l2.fsim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_fsim_30.XTick = [1e-3, 1e-2, 1e-1];
gca_fsim_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

ylim([0.9, 1])

xlabel('\lambda'); ylabel('FSIM - L (30% noise)')

subplot(1,2,2);
gca_nnz_30 = gca;
numel_30 = numel(O30.E);
semilogx(results_gl_l1.param1, results_gl_l1_err.nnz/numel_30, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2.param1, results_gl_l2_err.nnz/numel_30, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1.param1,    results_l1_err.nnz/numel_30,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2.param1,    results_l2_err.nnz/numel_30,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_nnz_30.XTick = [1e-3, 1e-2, 1e-1];
gca_nnz_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

ylim([0.28, 1])

xlabel('\lambda'); ylabel('Density - E (30% noise)')

%tightfig; 

fig_structure_30.Units = 'centimeters';
fig_structure_30.Position = [0 0 15.9 5.6];
fig_structure_30.PaperUnits = 'centimeters';
fig_structure_30.PaperPosition = [0 0 16 6];

saveas(fig_structure_30, 'synth_30_structure_mean', 'epsc');

%% FSIM and sparsity for 60% noise

fig_structure_60 = figure;

% subplot(1,2,1);
% gca_ssim_60 = gca;
% semilogx(results_gl_l1_60.param1, results_gl_l1_60.ssim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
% semilogx(results_gl_l2_60.param1, results_gl_l2_60.ssim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
% semilogx(results_l1_60.param1,    results_l1_60.ssim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
% semilogx(results_l2_60.param1,    results_l2_60.ssim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);
% 
% legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')
% 
% xlim([2e-4, 2e-1])
% gca_ssim_60.XTick = [1e-3, 1e-2, 1e-1];
% gca_ssim_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};
% 
% xlabel('\lambda'); ylabel('SSIM (60% noise)')

subplot(1,2,1);
gca_fsim_60 = gca;
semilogx(results_gl_l1_60.param1, results_gl_l1_60.fsim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2_60.param1, results_gl_l2_60.fsim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1_60.param1,    results_l1_60.fsim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2_60.param1,    results_l2_60.fsim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_fsim_60.XTick = [3e-3, 1e-2, 1e-1];
gca_fsim_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('FSIM - L (60% noise)')

subplot(1,2,2);
gca_nnz_60 = gca;
numel_60 = numel(O60.E);
semilogx(results_gl_l1.param1, results_gl_l1_60_err.nnz/numel_60, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2.param1, results_gl_l2_60_err.nnz/numel_60, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1.param1,    results_l1_60_err.nnz/numel_60,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2.param1,    results_l2_60_err.nnz/numel_60,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_nnz_60.XTick = [1e-3, 1e-2, 1e-1];
gca_nnz_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

ylim([0.58, 1])

xlabel('\lambda'); ylabel('Density - E (60% noise)')

%tightfig; 

fig_structure_60.Units = 'centimeters';
fig_structure_60.Position = [0 0 15.9 5.6];
fig_structure_60.PaperUnits = 'centimeters';
fig_structure_60.PaperPosition = [0 0 16 6];

saveas(fig_structure_60, 'synth_60_structure_mean', 'epsc');

%% Convergence plots

fig_cv = figure;

subplot(1,2,1);
gca_cv_30 = gca;
semilogy(best_gl_l1.rel_norm.iter, best_gl_l1.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_gl_l2.rel_norm.iter, best_gl_l2.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l1.rel_norm.iter, best_l1.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l2.rel_norm.iter, best_l2.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
gca_cv_30.ColorOrderIndex = 1;
semilogy(best_gl_l1.rel_norm.iter(1:5:end), best_gl_l1.rel_norm.err(1:5:end), 'v', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_gl_l2.rel_norm.iter(1:5:end), best_gl_l2.rel_norm.err(1:5:end), 'o', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l1.rel_norm.iter(1:5:end), best_l1.rel_norm.err(1:5:end), '+', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l2.rel_norm.iter(1:5:end), best_l2.rel_norm.err(1:5:end), 'x', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;

% legend plot
gca_cv_30.ColorOrderIndex = 1;
legend([legend_style('--',1.5,'v'), legend_style('--',1.5,'o'), ...
        legend_style('--',1.5,'+'), legend_style('--',1.5,'x')], ...
        {'GL L1', 'GL L2', 'Fro L1', 'Fro L2'}, ...
       'Location', 'northeast');

% legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

ylim([8e-8, 1.5])
gca_cv_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_cv_30.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('Number of iterations'); ylabel('Cv. criterion (30% noise)')

subplot(1,2,2);
gca_cv_60 = gca;
semilogy(best_gl_l1_60.rel_norm.iter, best_gl_l1_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_gl_l2_60.rel_norm.iter, best_gl_l2_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l1_60.rel_norm.iter, best_l1_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l2_60.rel_norm.iter, best_l2_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
gca_cv_60.ColorOrderIndex = 1;
semilogy(best_gl_l1_60.rel_norm.iter(1:5:end), best_gl_l1_60.rel_norm.err(1:5:end), 'v', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_gl_l2_60.rel_norm.iter(1:5:end), best_gl_l2_60.rel_norm.err(1:5:end), 'o', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l1_60.rel_norm.iter(1:5:end), best_l1_60.rel_norm.err(1:5:end), '+', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l2_60.rel_norm.iter(1:5:end), best_l2_60.rel_norm.err(1:5:end), 'x', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;

gca_cv_60.ColorOrderIndex = 1;
legend([legend_style('--',1.5,'v'), legend_style('--',1.5,'o'), ...
        legend_style('--',1.5,'+'), legend_style('--',1.5,'x')], ...
        {'GL L1', 'GL L2', 'Fro L1', 'Fro L2'}, ...
       'Location', 'northeast');

% legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

ylim([8e-8, 1.5])
gca_cv_60.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_cv_60.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('Number of iterations'); ylabel('Cv. criterion (60% noise)')

%tightfig; 

fig_cv.Units = 'centimeters';
fig_cv.Position = [0 0 15.9 5.6];
fig_cv.PaperUnits = 'centimeters';
fig_cv.PaperPosition = [0 0 16 6];

saveas(fig_cv, 'synth_cv_mean', 'epsc');

%% MSAM plots

fig_msam = figure;

subplot(1,2,1);
gca_msam_30 = gca;
loglog(results_gl_l1.param1, results_gl_l1.msam, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2.param1, results_gl_l2.msam, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1.param1, results_l1.msam, '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2.param1, results_l2.msam, '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5); hold on;

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_msam_30.XTick = [3e-3, 1e-2, 1e-1];

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')
gca_msam_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

%ylim([8e-8, 1.5])
%gca_msam_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_msam_30.YTickLabel = {'1e-8', '1e-6', '1e-4', '1e-2', '1'};

xlabel('\lambda'); ylabel('MSAM - L (30% noise)')

subplot(1,2,2);
gca_msam_60 = gca;
loglog(results_gl_l1_60.param1, results_gl_l1_60.msam, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2_60.param1, results_gl_l2_60.msam, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1_60.param1, results_l1_60.msam, '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2_60.param1, results_l2_60.msam, '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5); hold on;

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_msam_60.XTick = [3e-3, 1e-2, 1e-1];
gca_msam_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

%ylim([8e-8, 1.5])
%gca_msam_60.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_msam_60.YTickLabel = {'1e-8', '1e-6', '1e-4', '1e-2', '1'};

xlabel('\lambda'); ylabel('MSAM - L (60% noise)')

%tightfig; 

fig_msam.Units = 'centimeters';
fig_msam.Position = [0 0 15.9 5.6];
fig_msam.PaperUnits = 'centimeters';
fig_msam.PaperPosition = [0 0 16 6];

saveas(fig_msam, 'synth_msam_mean', 'epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_mean = figure;
subplot(1,2,1)
gca_mean_30 = gca;
loglog(results_gl_l1.param1, results_gl_l1_mean.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2.param1, results_gl_l2_mean.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1.param1,    results_l1_mean.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2.param1,    results_l2_mean.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_mean_30.XTick = [3e-3, 1e-2, 1e-1];
gca_mean_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('Rel. Err. - M (30% noise)')

ylim([1e-8, 1]);
gca_mean_30.YTick = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_mean_30.YTickLabel = {'1e-8','1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

subplot(1,2,2)
gca_mean_60 = gca;
loglog(results_gl_l1_60.param1, results_gl_l1_60_mean.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2_60.param1, results_gl_l2_60_mean.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1_60.param1,    results_l1_60_mean.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2_60.param1,    results_l2_60_mean.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_mean_60.XTick = [3e-3, 1e-2, 1e-1];
gca_mean_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('Rel. Err. - M (60% noise)')

ylim([4e-3, 1])
gca_mean_60.YTick = [1e-2, 1e-1, 1];
gca_mean_60.YTickLabel = {'1e-2', '1e-1', '1'};

%tightfig; 

fig_mean.Units = 'centimeters';
fig_mean.Position = [0 0 15.9 5.6];
fig_mean.PaperUnits = 'centimeters';
fig_mean.PaperPosition = [0 0 16 6];

saveas(fig_mean, 'synth_mean_rel_error', 'epsc');
