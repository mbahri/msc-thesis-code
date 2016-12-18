% [fg, tpr, fpr, auc, l_idx, r_idx] = select_best_foreground({D.L(:,:,50)}, O(:,:,50), true_mask, linspace(0, 1, 101));
% 
% elastic net: take weighted sum (convex) of l1 and l2 solutions?

%% Initialize and load the parameters
clear
close all


load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_rpca2d_gl_l1_sp_1.000000.mat')
[results_gl_l1, best_gl_l1] = distilate(results_mehdi.rpca2d_gl_l1, 'lambda');
[results_gl_l1_err, best_gl_l1_err] = distilate_err(results_mehdi.rpca2d_gl_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_rpca2d_gl_l2_sp_1.000000.mat')
[results_gl_l2, best_gl_l2] = distilate(results_mehdi.rpca2d_gl_l2, 'lambda');
[results_gl_l2_err, best_gl_l2_err] = distilate_err(results_mehdi.rpca2d_gl_l2, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_rpca2d_l1_sp_1.000000.mat')
[results_l1, best_l1] = distilate(results_mehdi.rpca2d_l1, 'lambda');
[results_l1_err, best_l1_err] = distilate_err(results_mehdi.rpca2d_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_30_mehdi_rpca2d_l2_sp_1.000000.mat')
[results_l2, best_l2] = distilate(results_mehdi.rpca2d_l2, 'lambda');
[results_l2_err, best_l2_err] = distilate_err(results_mehdi.rpca2d_l2, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_rpca2d_gl_l1_sp_1.000000.mat')
[results_gl_l1_60, best_gl_l1_60] = distilate(results_mehdi.rpca2d_gl_l1, 'lambda');
[results_gl_l1_60_err, best_gl_l1_60_err] = distilate_err(results_mehdi.rpca2d_gl_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_rpca2d_gl_l2_sp_1.000000.mat')
[results_gl_l2_60, best_gl_l2_60] = distilate(results_mehdi.rpca2d_gl_l2, 'lambda');
[results_gl_l2_60_err, best_gl_l2_60_err] = distilate_err(results_mehdi.rpca2d_gl_l2, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_rpca2d_l1_sp_1.000000.mat')
[results_l1_60, best_l1_60] = distilate(results_mehdi.rpca2d_l1, 'lambda');
[results_l1_60_err, best_l1_60_err] = distilate_err(results_mehdi.rpca2d_l1, 'lambda');

load('/vol/bitbucket/mb2215/Thesis/Results/Synthetic/results_synthetic_60_mehdi_rpca2d_l2_sp_1.000000.mat')
[results_l2_60, best_l2_60] = distilate(results_mehdi.rpca2d_l2, 'lambda');
[results_l2_60_err, best_l2_60_err] = distilate_err(results_mehdi.rpca2d_l2, 'lambda');


%% Plot PSNR and relative error for 30% noise
fig_err_psnr_30 = figure

subplot(1,2,1);
gca_psnr_30 = gca;
semilogx(results_gl_l1.param1, results_gl_l1.psnr, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2.param1, results_gl_l2.psnr, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1.param1,    results_l1.psnr,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2.param1,    results_l2.psnr,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_psnr_30.XTick = [1e-3, 1e-2, 1e-1];
gca_psnr_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('PSNR (30% noise)')

subplot(1,2,2)
gca_err_30 = gca;
loglog(results_gl_l1.param1, results_gl_l1.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2.param1, results_gl_l2.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1.param1,    results_l1.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2.param1,    results_l2.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_err_30.XTick = [1e-3, 1e-2, 1e-1];
gca_err_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 30%
ylim([5e-8, 1])
gca_err_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_err_30.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('\lambda'); ylabel('Relative error (30% noise)')

tightfig;

fig_err_psnr_30.Units = 'centimeters';
fig_err_psnr_30.Position = [0 0 15.9 8.8];
fig_err_psnr_30.PaperUnits = 'centimeters';
fig_err_psnr_30.PaperPosition = [0 0 16 9];

saveas(fig_err_psnr_30, 'synth_30_err_psnr', 'epsc');

%% Plot PSNR and relative error for 60% noise
fig_err_psnr_60 = figure

subplot(1,2,1);
gca_psnr_60 = gca;
semilogx(results_gl_l1_60.param1, results_gl_l1_60.psnr, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2_60.param1, results_gl_l2_60.psnr, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1_60.param1,    results_l1_60.psnr,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2_60.param1,    results_l2_60.psnr,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_psnr_60.XTick = [1e-3, 1e-2, 1e-1];
gca_psnr_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 60%
ylim([18, 160]);

xlabel('\lambda'); ylabel('PSNR (60% noise)')

subplot(1,2,2)
gca_err_60 = gca;
loglog(results_gl_l1_60.param1, results_gl_l1_60.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2_60.param1, results_gl_l2_60.rel_norm, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1_60.param1,    results_l1_60.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2_60.param1,    results_l2_60.rel_norm,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_err_60.XTick = [1e-3, 1e-2, 1e-1];
gca_err_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

% 60%
ylim([1e-7, 1]);
gca_err_60.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_err_60.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('\lambda'); ylabel('Relative error (60% noise)')

tightfig;

fig_err_psnr_60.Units = 'centimeters';
fig_err_psnr_60.Position = [0 0 15.9 8.8];
fig_err_psnr_60.PaperUnits = 'centimeters';
fig_err_psnr_60.PaperPosition = [0 0 16 9];

saveas(fig_err_psnr_60, 'synth_60_err_psnr', 'epsc');

%% SSIM and FSIM for 30% noise

fig_structure_30 = figure

subplot(1,2,1);
gca_ssim_30 = gca;
semilogx(results_gl_l1.param1, results_gl_l1.ssim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2.param1, results_gl_l2.ssim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1.param1,    results_l1.ssim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2.param1,    results_l2.ssim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_ssim_30.XTick = [1e-3, 1e-2, 1e-1];
gca_ssim_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('SSIM (30% noise)')

subplot(1,2,2);
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

xlabel('\lambda'); ylabel('FSIM (30% noise)')

tightfig;

fig_structure_30.Units = 'centimeters';
fig_structure_30.Position = [0 0 15.9 8.8];
fig_structure_30.PaperUnits = 'centimeters';
fig_structure_30.PaperPosition = [0 0 16 9];

saveas(fig_structure_30, 'synth_30_ssim_fsim', 'epsc');

%% SSIM and FSIM for 60% noise

fig_structure_60 = figure

subplot(1,2,1);
gca_ssim_60 = gca;
semilogx(results_gl_l1_60.param1, results_gl_l1_60.ssim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2_60.param1, results_gl_l2_60.ssim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1_60.param1,    results_l1_60.ssim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2_60.param1,    results_l2_60.ssim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_ssim_60.XTick = [1e-3, 1e-2, 1e-1];
gca_ssim_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('SSIM (60% noise)')

subplot(1,2,2);
gca_fsim_60 = gca;
semilogx(results_gl_l1_60.param1, results_gl_l1_60.fsim, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
semilogx(results_gl_l2_60.param1, results_gl_l2_60.fsim, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
semilogx(results_l1_60.param1,    results_l1_60.fsim,    '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
semilogx(results_l2_60.param1,    results_l2_60.fsim,    '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5);

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'southeast')

xlim([2e-4, 2e-1])
gca_fsim_60.XTick = [3e-3, 1e-2, 1e-1];
gca_fsim_60.XTickLabel = {'1e-3', '1e-2', '1e-1'};

xlabel('\lambda'); ylabel('FSIM (60% noise)')

tightfig;

fig_structure_60.Units = 'centimeters';
fig_structure_60.Position = [0 0 15.9 8.8];
fig_structure_60.PaperUnits = 'centimeters';
fig_structure_60.PaperPosition = [0 0 16 9];

saveas(fig_structure_60, 'synth_60_ssim_fsim', 'epsc');

%% Convergence plots

fig_cv = figure

subplot(1,2,1);
gca_cv_30 = gca
semilogy(best_gl_l1.rel_norm.iter, best_gl_l1.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_gl_l2.rel_norm.iter, best_gl_l2.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l1.rel_norm.iter, best_l1.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l2.rel_norm.iter, best_l2.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
gca_cv_30.ColorOrderIndex = 1;
semilogy(best_gl_l1.rel_norm.iter(1:5:end), best_gl_l1.rel_norm.err(1:5:end), 'v', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_gl_l2.rel_norm.iter(1:5:end), best_gl_l2.rel_norm.err(1:5:end), 'o', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l1.rel_norm.iter(1:5:end), best_l1.rel_norm.err(1:5:end), '+', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l2.rel_norm.iter(1:5:end), best_l2.rel_norm.err(1:5:end), 'x', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

ylim([8e-8, 1.5])
gca_cv_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_cv_30.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('Number of iterations'); ylabel('Max. rel. error (30% noise)')

subplot(1,2,2);
gca_cv_60 = gca
semilogy(best_gl_l1_60.rel_norm.iter, best_gl_l1_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_gl_l2_60.rel_norm.iter, best_gl_l2_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l1_60.rel_norm.iter, best_l1_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
semilogy(best_l2_60.rel_norm.iter, best_l2_60.rel_norm.err, '--', 'LineWidth', 1.5); hold on;
gca_cv_60.ColorOrderIndex = 1;
semilogy(best_gl_l1_60.rel_norm.iter(1:5:end), best_gl_l1_60.rel_norm.err(1:5:end), 'v', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_gl_l2_60.rel_norm.iter(1:5:end), best_gl_l2_60.rel_norm.err(1:5:end), 'o', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l1_60.rel_norm.iter(1:5:end), best_l1_60.rel_norm.err(1:5:end), '+', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogy(best_l2_60.rel_norm.iter(1:5:end), best_l2_60.rel_norm.err(1:5:end), 'x', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

ylim([8e-8, 1.5])
gca_cv_60.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_cv_60.YTickLabel = {'1e-7', '1e-6', '1e-5', '1e-5', '1e-3', '1e-2', '1e-1', '1'};

xlabel('Number of iterations'); ylabel('Max. rel. error (60% noise)')

tightfig;

fig_cv.Units = 'centimeters';
fig_cv.Position = [0 0 15.9 8.6];
fig_cv.PaperUnits = 'centimeters';
fig_cv.PaperPosition = [0 0 16 9];

saveas(fig_cv, 'synth_cv', 'epsc');

%% MSAM plots

fig_msam = figure

subplot(1,2,1);
gca_msam_30 = gca
loglog(results_gl_l1.param1, results_gl_l1.msam, '--', 'LineWidth', 1.5, 'Marker', 'v', 'MarkerSize', 5); hold on;
loglog(results_gl_l2.param1, results_gl_l2.msam, '--', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5); hold on;
loglog(results_l1.param1, results_l1.msam, '--', 'LineWidth', 1.5, 'Marker', '+', 'MarkerSize', 5); hold on;
loglog(results_l2.param1, results_l2.msam, '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 5); hold on;

legend('GL L1', 'GL L2', 'Fro L1', 'Fro L2', 'Location', 'northeast')

xlim([2e-4, 2e-1])
gca_msam_30.XTick = [3e-3, 1e-2, 1e-1];
gca_msam_30.XTickLabel = {'1e-3', '1e-2', '1e-1'};

%ylim([8e-8, 1.5])
%gca_msam_30.YTick = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
gca_msam_30.YTickLabel = {'1e-8', '1e-6', '1e-4', '1e-2', '1'};

xlabel('\lambda'); ylabel('MSAM (30% noise)')

subplot(1,2,2);
gca_msam_60 = gca
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

xlabel('\lambda'); ylabel('MSAM (60% noise)')

tightfig;

fig_msam.Units = 'centimeters';
fig_msam.Position = [0 0 15.9 8.6];
fig_msam.PaperUnits = 'centimeters';
fig_msam.PaperPosition = [0 0 16 9];

saveas(fig_msam, 'synth_msam', 'epsc');
