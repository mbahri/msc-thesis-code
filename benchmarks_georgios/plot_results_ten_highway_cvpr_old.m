% Plots the results of the background subtraction experiment on the highway
% video using tensor methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
close all;
% clc;

% choose data
dims = [48 64 400];
method_type = 'tensors';

% load data
path = '/vol/bitbucket/mb2215/experiments/highway';
file1 = sprintf('%dx%dx%d_no_noise.mat', dims);
file2 = sprintf('%dx%dx%d_no_noise_part2.mat', dims);
file3 = sprintf('%dx%dx%d_no_noise_cauchy.mat', dims);
file4 = sprintf('%dx%dx%d_no_noise_welsh.mat', dims);
file5 = sprintf('%dx%dx%d_no_noise_nctrpca.mat', dims);

load(fullfile(path, method_type, file1));
load(fullfile(path, method_type, file2));
load(fullfile(path, method_type, file3));
load(fullfile(path, method_type, file4));
load(fullfile(path, method_type, file5));

load(fullfile(path, 'params_no_noise.mat'));
load(fullfile(path, 'params_no_noise_part2.mat'));
load(fullfile(path, 'params_no_noise_cauchy.mat'));
load(fullfile(path, 'params_no_noise_welsh.mat'));
load(fullfile(path, 'params_no_noise_nctrpca.mat'));

results.cauchy_st = results_cauchy.cauchy_st;
results.welsh_st = results_welsh.welsh_st;
results.brtf = results_part2.brtf;
results.horpca_s = results_part2.horpca_s;
results.nctrpca = results_nctrpca.nctrpca;


% plot parameters
axis_fontsize = 14;
linewidth = 2;
markersize_big = 5;
markersize_sml = 2;
colour = distinguishable_colors(20);
% colour = colour([1 5 4 3 7 6 2 8], :);
colour = colour(1:20, :);
marker = {'d', 'x', 's', '>', '<', '+', '*', 'p'};
linestyle = ':';
position = [0 0 26 11];
units = 'centimeters';
paperposition = [0 0 26.5 11.5];

% display pics params
% load('../data/highway.mat');
load_highway
rowidx = [1 2];
colidx = 3;

imsize = dims(rowidx);
range = [min(X(:)) max(X(:))];

% figure save path
savepath = fullfile('/vol/bitbucket/mb2215/Thesis', 'results/cvpr', 'highway');

% % squeeze matrices
% results.brpca.err  = squeeze(results.brpca.err );
% results.brpca.time = squeeze(results.brpca.time);
% results.brpca.info = squeeze(results.brpca.info);
% results.brpca.A    = squeeze(results.brpca.A   );
% 
% results.irpca.sub.err  = squeeze(results.irpca.sub.err );
% results.irpca.sub.time = squeeze(results.irpca.sub.time);
% results.irpca.sub.info = squeeze(results.irpca.sub.info);
% results.irpca.sub.A    = squeeze(results.irpca.sub.A   );
% 
% results.irpca.lin.err  = squeeze(results.irpca.lin.err );
% results.irpca.lin.time = squeeze(results.irpca.lin.time);
% results.irpca.lin.info = squeeze(results.irpca.lin.info);
% results.irpca.lin.A    = squeeze(results.irpca.lin.A   );
% 
% results.orpca.err  = squeeze(results.orpca.err );
% results.orpca.time = squeeze(results.orpca.time);
% results.orpca.info = squeeze(results.orpca.info);
% results.orpca.A    = squeeze(results.orpca.A   );
% 
% results.rcpd.sub.err  = squeeze(results.rcpd.sub.err );
% results.rcpd.sub.time = squeeze(results.rcpd.sub.time);
% results.rcpd.sub.info = squeeze(results.rcpd.sub.info);
% results.rcpd.sub.A    = squeeze(results.rcpd.sub.A   );
% 
% results.rcpd.lin.err  = squeeze(results.rcpd.lin.err );
% results.rcpd.lin.time = squeeze(results.rcpd.lin.time);
% results.rcpd.lin.info = squeeze(results.rcpd.lin.info);
% results.rcpd.lin.A    = squeeze(results.rcpd.lin.A   );

%% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_tpr   = containers.Map;
best_fpr   = containers.Map;
best_auc   = containers.Map;
best_foreground = containers.Map;

true_mask = GT > min(GT(:)) + 0.1;
threshold = linspace(0, 1, 101);

% [best_foreground('rpca'     ), best_tpr('rpca'     ), best_fpr('rpca'     ), best_auc('rpca'     ), best_l_idx('rpca'     )                         ] = select_best_foreground(results. rpca.    A, X, true_mask, threshold);
% [best_foreground('brpca'    ), best_tpr('brpca'    ), best_fpr('brpca'    ), best_auc('brpca'    ), best_l_idx('brpca'    ), best_r_idx('brpca'    )] = select_best_foreground(results.brpca.    A, X, true_mask, threshold);
% [best_foreground('irpca_sub'), best_tpr('irpca_sub'), best_fpr('irpca_sub'), best_auc('irpca_sub'), best_l_idx('irpca_sub'), best_r_idx('irpca_sub')] = select_best_foreground(results.irpca.sub.A, X, true_mask, threshold);
% [best_foreground('irpca_lin'), best_tpr('irpca_lin'), best_fpr('irpca_lin'), best_auc('irpca_lin'), best_l_idx('irpca_lin'), best_r_idx('irpca_lin')] = select_best_foreground(results.irpca.lin.A, X, true_mask, threshold);
% [best_foreground('orpca'    ), best_tpr('orpca'    ), best_fpr('orpca'    ), best_auc('orpca'    ), best_l_idx('orpca'    ), best_r_idx('orpca'    )] = select_best_foreground(results.orpca.    A, X, true_mask, threshold);
% [best_foreground('rcpd_sub' ), best_tpr('rcpd_sub' ), best_fpr('rcpd_sub' ), best_auc('rcpd_sub' ), best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )] = select_best_foreground(results. rcpd.sub.A, X, true_mask, threshold);
% [best_foreground('rcpd_lin' ), best_tpr('rcpd_lin' ), best_fpr('rcpd_lin' ), best_auc('rcpd_lin' ), best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )] = select_best_foreground(results. rcpd.lin.A, X, true_mask, threshold);
% 
% [best_foreground('rpca2d_l2' ), best_tpr('rpca2d_l2' ), best_fpr('rpca2d_l2' ), best_auc('rpca2d_l2' ), best_l_idx('rpca2d_l2' ), best_r_idx('rpca2d_l2' )] = select_best_foreground(results.rpca2d_l2.A, X, true_mask, threshold);
[best_foreground('rpca2d_l1' ), best_tpr('rpca2d_l1' ), best_fpr('rpca2d_l1' ), best_auc('rpca2d_l1' ), best_l_idx('rpca2d_l1' ), best_r_idx('rpca2d_l1' )] = select_best_foreground(results.rpca2d_l1.A, X, true_mask, threshold);
% [best_foreground('rpca2d_gl_l2' ), best_tpr('rpca2d_gl_l2' ), best_fpr('rpca2d_gl_l2' ), best_auc('rpca2d_gl_l2' ), best_l_idx('rpca2d_gl_l2' ), best_r_idx('rpca2d_gl_l2' )] = select_best_foreground(results.rpca2d_gl_l2.A, X, true_mask, threshold);
% [best_foreground('rpca2d_gl_l1' ), best_tpr('rpca2d_gl_l1' ), best_fpr('rpca2d_gl_l1' ), best_auc('rpca2d_gl_l1' ), best_l_idx('rpca2d_gl_l1' ), best_r_idx('rpca2d_gl_l1' )] = select_best_foreground(results.rpca2d_gl_l1.A, X, true_mask, threshold);

[best_foreground('tnn' ), best_tpr('tnn' ), best_fpr('tnn' ), best_auc('tnn' ), best_l_idx('tnn' ), best_r_idx('tnn' )] = select_best_foreground(results.tnn.A, X, true_mask, threshold);
[best_foreground('tsvd' ), best_tpr('tsvd' ), best_fpr('tsvd' ), best_auc('tsvd' ), best_l_idx('tsvd' ), best_r_idx('tsvd' )] = select_best_foreground(results.tsvd.A, X, true_mask, threshold);

[best_foreground('cauchy_st' ), best_tpr('cauchy_st' ), best_fpr('cauchy_st' ), best_auc('cauchy_st' ), best_l_idx('cauchy_st' ), best_r_idx('cauchy_st' )] = select_best_foreground(results.cauchy_st.A, X, true_mask, threshold);
[best_foreground('welsh_st' ), best_tpr('welsh_st' ), best_fpr('welsh_st' ), best_auc('welsh_st' ), best_l_idx('welsh_st' ), best_r_idx('welsh_st' )] = select_best_foreground(results.welsh_st.A, X, true_mask, threshold);
[best_foreground('horpca_s' ), best_tpr('horpca_s' ), best_fpr('horpca_s' ), best_auc('horpca_s' ), best_l_idx('horpca_s' ), best_r_idx('horpca_s' )] = select_best_foreground(results.horpca_s.A, X, true_mask, threshold);
% [best_foreground('brtf' ), best_tpr('brtf' ), best_fpr('brtf' ), best_auc('brtf' ), best_l_idx('brtf' ), best_r_idx('brtf' )] = select_best_foreground(results.brtf.A, X, true_mask, threshold);

[best_foreground('nctrpca' ), best_tpr('nctrpca' ), best_fpr('nctrpca' ), best_auc('nctrpca' ), best_l_idx('nctrpca' ), best_r_idx('nctrpca' )] = select_best_foreground(results.nctrpca.A, X, true_mask, threshold);

maxiter = max([results. rpca2d_l1.info{best_l_idx('rpca2d_l1' ), best_r_idx('rpca2d_l1' )}.iter(end) ...
               results. tnn.info{best_l_idx('tnn' ), best_r_idx('tnn' )}.iter(end) ...
               results. horpca_s.info{best_l_idx('horpca_s' ), best_r_idx('horpca_s' )}.iter(end) ...
               results. nctrpca.info{best_l_idx('nctrpca' ), best_r_idx('nctrpca' )}.iter(end) ...
               ] );


%% ROC curves
cnt = 1;

fig = figure('Position', position); 
% fpr = norminv(best_fpr('rpca'     )); fnr = norminv(1 - best_tpr('rpca'     )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
% fpr = norminv(best_fpr('brpca'    )); fnr = norminv(1 - best_tpr('brpca'    )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% % fpr = norminv(best_fpr('irpca_sub')); fnr = norminv(1 - best_tpr('irpca_sub')); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% % fpr = norminv(best_fpr('irpca_lin')); fnr = norminv(1 - best_tpr('irpca_lin')); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% fpr = norminv(best_fpr('orpca'    )); fnr = norminv(1 - best_tpr('orpca'    )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% fpr = norminv(best_fpr('rcpd_sub' )); fnr = norminv(1 - best_tpr('rcpd_sub' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% fpr = norminv(best_fpr('rcpd_lin' )); fnr = norminv(1 - best_tpr('rcpd_lin' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% 
% fpr = norminv(best_fpr('rpca2d_l2' )); fnr = norminv(1 - best_tpr('rpca2d_l2' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{2}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('rpca2d_l1' )); fnr = norminv(1 - best_tpr('rpca2d_l1' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{2}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
% fpr = norminv(best_fpr('rpca2d_gl_l2' )); fnr = norminv(1 - best_tpr('rpca2d_gl_l2' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{2}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% fpr = norminv(best_fpr('rpca2d_gl_l1' )); fnr = norminv(1 - best_tpr('rpca2d_gl_l1' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{2}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('tnn' )); fnr = norminv(1 - best_tpr('tnn' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{3}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('tsvd' )); fnr = norminv(1 - best_tpr('tsvd' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{4}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('horpca_s' )); fnr = norminv(1 - best_tpr('horpca_s' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{5}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
% fpr = norminv(best_fpr('brtf' )); fnr = norminv(1 - best_tpr('brtf' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{6}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('cauchy_st' )); fnr = norminv(1 - best_tpr('cauchy_st' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{7}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('welsh_st' )); fnr = norminv(1 - best_tpr('welsh_st' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{7}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('nctrpca' )); fnr = norminv(1 - best_tpr('nctrpca' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{8}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth);

ylabel('False Negative Rate', 'FontSize', axis_fontsize);
xlabel('False Positive Rate', 'FontSize', axis_fontsize);
legend('KDRSDL (Proposed)', 'Tensor RPCA (CPVR2016)', 'Tensor RPCA (CPVR2014)', ...
    'Tensor HORPCA-S (Goldfarb & Qin)', 'Cauchy Soft-Thresholding', 'Welsh Soft-Thresholding', 'NC TRPCA', ...
    'Location', 'NorthEastOutside');

set(gca, 'FontSize', axis_fontsize);
ppxtick = @(x) ( strrep(sprintf('%.0e',x), 'e-0', 'e-') );
xtick = [10.^(-4:-1) 0.5]; set(gca, 'XTick', norminv(xtick)); set(gca, 'XTickLabel', arrayfun(ppxtick, xtick, 'UniformOutput', 0));
ytick = 0.1 : 0.2 : 0.9;   set(gca, 'YTick', norminv(ytick)); set(gca, 'YTickLabel', {ytick});
xlim(norminv([0.0000007 0.76])); 
ylim(norminv([0.05 0.992]));
set(fig, 'PaperPositionMode', 'auto'); 
saveas(fig, fullfile(savepath, 'ten_highway_roc_curves.eps'), 'epsc');

fig.PaperOrientation = 'landscape';
fig.Renderer = 'painters';
fig.Units = units;
fig.Position = position;
fig.PaperUnits = units;
fig.PaperPosition = paperposition;
print(fig, '-depsc', fullfile(savepath, 'ten_highway_roc_curves'));

%% best results
% fprintf('Tensor RPCA (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('rpca')));
% fprintf('\t time = %g sec \n', results.rpca.time(best_l_idx('rpca')));
% fprintf('\t auc = %g \n', best_auc('rpca'));
% fprintf('\t iterations = %d \n\n', results.rpca.info{best_l_idx('rpca')}.iter(end));
% 
% fprintf('Tensor BRPCA (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('brpca')));
% fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('brpca')));
% fprintf('\t time = %g sec \n', results.brpca.time(best_l_idx('brpca'),best_r_idx('brpca')));
% fprintf('\t auc = %g \n', best_auc('brpca'));
% fprintf('\t iterations = %d \n\n', results.brpca.info{best_l_idx('brpca'),best_r_idx('brpca')}.iter(end));
% 
% fprintf('Tensor IRPCA (sub) (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_sub')));
% fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('irpca_sub')));
% fprintf('\t time = %g sec \n', results.irpca.sub.time(best_l_idx('irpca_sub'),best_r_idx('irpca_sub')));
% fprintf('\t auc = %g \n', best_auc('irpca_sub'));
% fprintf('\t iterations = %d \n\n', results.irpca.sub.info{best_l_idx('irpca_sub'),best_r_idx('irpca_sub')}.iter(end));
% 
% fprintf('Tensor IRPCA (lin) (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('irpca_lin')));
% fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('irpca_lin')));
% fprintf('\t time = %g sec \n', results.irpca.lin.time(best_l_idx('irpca_lin'),best_r_idx('irpca_lin')));
% fprintf('\t auc = %g \n', best_auc('irpca_lin'));
% fprintf('\t iterations = %d \n\n', results.irpca.lin.info{best_l_idx('irpca_lin'),best_r_idx('irpca_lin')}.iter(end));
% 
% fprintf('Tensor RHOSVD (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('orpca')));
% fprintf('\t nrank [%%] = %g \n', ten_nrank_perc(best_r_idx('orpca')));
% fprintf('\t time = %g sec \n', results.orpca.time(best_l_idx('orpca'),best_r_idx('orpca')));
% fprintf('\t auc = %g \n', best_auc('orpca'));
% fprintf('\t iterations = %d \n\n', results.orpca.info{best_l_idx('orpca'),best_r_idx('orpca')}.iter(end));
% 
% fprintf('Tensor RCPD (sub) (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('rcpd_sub')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('rcpd_sub')));
% fprintf('\t time = %g sec \n', results.rcpd.sub.time(best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')));
% fprintf('\t auc = %g \n', best_auc('rcpd_sub'));
% fprintf('\t iterations = %d \n\n', results.rcpd.sub.info{best_l_idx('rcpd_sub'),best_r_idx('rcpd_sub')}.iter(end));
% 
% fprintf('Tensor RCPD (lin) (Papamakarios et. al.): \n');
% fprintf('\t lambda = %g \n', lambda(best_l_idx('rcpd_lin')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('rcpd_lin')));
% fprintf('\t time = %g sec \n', results.rcpd.lin.time(best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')));
% fprintf('\t auc = %g \n', best_auc('rcpd_lin'));
% fprintf('\t iterations = %d \n\n', results.rcpd.lin.info{best_l_idx('rcpd_lin'),best_r_idx('rcpd_lin')}.iter(end));
% 
% fprintf('RPCA2D L2 (Bahri et. al.): \n');
% fprintf('\t lambda = %g \n', params_2d.lambda(best_l_idx('rpca2d_l2')));
% % fprintf('\t rank = %g \n', ten_rank(best_r_idx('rpca2d_l2')));
% fprintf('\t time = %g sec \n', results.rpca2d_l2.time(best_l_idx('rpca2d_l2'),best_r_idx('rpca2d_l2')));
% fprintf('\t auc = %g \n', best_auc('rpca2d_l2'));
% fprintf('\t iterations = %d \n\n', results.rpca2d_l2.info{best_l_idx('rpca2d_l2'),best_r_idx('rpca2d_l2')}.iter(end));

fprintf('RPCA2D L1 (Bahri et. al.): \n');
fprintf('\t lambda = %g \n', params_2d.lambda(best_l_idx('rpca2d_l1')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('rpca2d_l1')));
fprintf('\t time = %g sec \n', results.rpca2d_l2.time(best_l_idx('rpca2d_l1'),best_r_idx('rpca2d_l1')));
fprintf('\t auc = %g \n', best_auc('rpca2d_l1'));
fprintf('\t iterations = %d \n\n', results.rpca2d_l1.info{best_l_idx('rpca2d_l1'),best_r_idx('rpca2d_l1')}.iter(end));

% fprintf('RPCA2D GL L2 (Bahri et. al.): \n');
% fprintf('\t lambda = %g \n', params_2d_gl.lambda(best_l_idx('rpca2d_gl_l2')));
% fprintf('\t alpha = %g \n', params_2d_gl.p(best_r_idx('rpca2d_gl_l2')));
% fprintf('\t time = %g sec \n', results.rpca2d_gl_l2.time(best_l_idx('rpca2d_gl_l2'),best_r_idx('rpca2d_gl_l2')));
% fprintf('\t auc = %g \n', best_auc('rpca2d_gl_l2'));
% fprintf('\t iterations = %d \n\n', results.rpca2d_gl_l2.info{best_l_idx('rpca2d_gl_l2'),best_r_idx('rpca2d_gl_l2')}.iter(end));
% 
% fprintf('RPCA2D GL L1 (Bahri et. al.): \n');
% fprintf('\t lambda = %g \n', params_2d_gl.lambda(best_l_idx('rpca2d_gl_l1')));
% fprintf('\t alpha = %g \n', params_2d_gl.p(best_r_idx('rpca2d_gl_l1')));
% fprintf('\t time = %g sec \n', results.rpca2d_gl_l1.time(best_l_idx('rpca2d_gl_l1'),best_r_idx('rpca2d_gl_l1')));
% fprintf('\t auc = %g \n', best_auc('rpca2d_gl_l1'));
% fprintf('\t iterations = %d \n\n', results.rpca2d_gl_l1.info{best_l_idx('rpca2d_gl_l1'),best_r_idx('rpca2d_gl_l1')}.iter(end));

fprintf('Tensor RPCA (CPVR2016): \n');
fprintf('\t lambda = %g \n', params_tnn.lambda(best_l_idx('tnn')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('tnn')));
fprintf('\t time = %g sec \n', results.tnn.time(best_l_idx('tnn'),best_r_idx('tnn')));
fprintf('\t auc = %g \n', best_auc('tnn'));
fprintf('\t iterations = %d \n\n', results.tnn.info{best_l_idx('tnn'),best_r_idx('tnn')}.iter(end));

fprintf('Tensor RPCA (CPVR2014): \n');
fprintf('\t lambda = %g \n', params_tsvd.lambda(best_l_idx('tsvd')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('tsvd')));
fprintf('\t time = %g sec \n', results.tsvd.time(best_l_idx('tsvd'),best_r_idx('tsvd')));
fprintf('\t auc = %g \n', best_auc('tsvd'));
fprintf('\t iterations = %d \n\n', results.tsvd.info{best_l_idx('tsvd'),best_r_idx('tsvd')}.iter(end));

fprintf('Tensor HORPCA-S (Goldfarb & Qin): \n');
fprintf('\t lambda = %g \n', params_horpca.lambda(best_l_idx('horpca_s')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('horpca_s')));
fprintf('\t time = %g sec \n', results.horpca_s.time(best_l_idx('horpca_s'),best_r_idx('horpca_s')));
fprintf('\t auc = %g \n', best_auc('horpca_s'));
fprintf('\t iterations = %d \n\n', results.horpca_s.info{best_l_idx('horpca_s'),best_r_idx('horpca_s')}.iter(end));

% fprintf('BRTF: \n');
% fprintf('\t initVar = %g \n', params_brtf.lambda(best_l_idx('brtf')));
% fprintf('\t maxRank = %g \n', 100);
% fprintf('\t time = %g sec \n', results.brtf.time(best_l_idx('brtf'),best_r_idx('brtf')));
% fprintf('\t auc = %g \n', best_auc('brtf'));
% fprintf('\t iterations = %d \n\n', results.brtf.info{best_l_idx('brtf'),best_r_idx('brtf')}.iter(end));

fprintf('Cauchy Soft-Thresholding: \n');
fprintf('\t sigma = %g \n', params_mest.p(best_r_idx('cauchy_st')));
fprintf('\t alpha = %g \n', params_mest.lambda(best_l_idx('cauchy_st')));
% fprintf('\t rank = %g \n', 0);
fprintf('\t time = %g sec \n', results.cauchy_st.time(best_l_idx('cauchy_st'),best_r_idx('cauchy_st')));
fprintf('\t auc = %g \n', best_auc('cauchy_st'));
fprintf('\t iterations = %d \n\n', results.cauchy_st.info{best_l_idx('cauchy_st'),best_r_idx('cauchy_st')}.iter(end));

fprintf('Welsh Soft-Thresholding: \n');
fprintf('\t sigma = %g \n', params_mest.p(best_r_idx('welsh_st')));
fprintf('\t alpha = %g \n', params_mest.lambda(best_l_idx('welsh_st')));
% fprintf('\t rank = %g \n', 0);
fprintf('\t time = %g sec \n', results.welsh_st.time(best_l_idx('welsh_st'),best_r_idx('welsh_st')));
fprintf('\t auc = %g \n', best_auc('welsh_st'));
fprintf('\t iterations = %d \n\n', results.welsh_st.info{best_l_idx('welsh_st'),best_r_idx('welsh_st')}.iter(end));

fprintf('NC TRPCA: \n');
fprintf('\t threshold = %g \n', params_nctrpca.p(best_r_idx('nctrpca')));
fprintf('\t rank = %g \n', params_nctrpca.lambda(best_l_idx('nctrpca')));
fprintf('\t time = %g sec \n', results.nctrpca.time(best_l_idx('nctrpca'),best_r_idx('nctrpca')));
fprintf('\t auc = %g \n', best_auc('nctrpca'));
fprintf('\t iterations = %d \n\n', results.nctrpca.info{best_l_idx('nctrpca'),best_r_idx('nctrpca')}.iter(end));

%% save images
% im_idx = [130 164 315 355];
% imformat = 'png';
% 
% for i = 1:length(im_idx)
%     imfile = sprintf('ten_highway_original_%d.%s', im_idx(i), imformat); im = X        (:,:,im_idx(i)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_highway_mask_%d.%s'    , im_idx(i), imformat); im = true_mask(:,:,im_idx(i)); imwrite(im, fullfile(savepath, imfile), imformat);
%     
%     imfile = sprintf('ten_bg_rpca_%d.%s'     , im_idx(i), imformat); im = results. rpca.    A{best_l_idx('rpca'     )                         }(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_bg_brpca_%d.%s'    , im_idx(i), imformat); im = results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca'    )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_bg_irpca_sub_%d.%s', im_idx(i), imformat); im = results.irpca.sub.A{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_bg_irpca_lin_%d.%s', im_idx(i), imformat); im = results.irpca.lin.A{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_bg_orpca_%d.%s'    , im_idx(i), imformat); im = results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca'    )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_bg_rcpd_sub_%d.%s' , im_idx(i), imformat); im = results. rcpd.sub.A{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_bg_rcpd_lin_%d.%s' , im_idx(i), imformat); im = results. rcpd.lin.A{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}(:,:,im_idx(i)); imwrite(double(im), fullfile(savepath, imfile), imformat);
%     
%     imfile = sprintf('ten_fg_rpca_%d.%s'     , im_idx(i), imformat); im = best_foreground('rpca'     ); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_fg_brpca_%d.%s'    , im_idx(i), imformat); im = best_foreground('brpca'    ); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_fg_irpca_sub_%d.%s', im_idx(i), imformat); im = best_foreground('irpca_sub'); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_fg_irpca_lin_%d.%s', im_idx(i), imformat); im = best_foreground('irpca_lin'); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_fg_orpca_%d.%s'    , im_idx(i), imformat); im = best_foreground('orpca'    ); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_fg_rcpd_sub_%d.%s' , im_idx(i), imformat); im = best_foreground('rcpd_sub' ); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_fg_rcpd_lin_%d.%s' , im_idx(i), imformat); im = best_foreground('rcpd_lin' ); im = im(:,:,im_idx(i)); im = im - min(im(:)); im = im / max(im(:)); imwrite(im, fullfile(savepath, imfile), imformat);
% end
% 
% %% save movies
% movformat = 'avi';
% movfile = sprintf('ten_highway.%s', movformat); arr2movie(cat(2, X, true_mask), fullfile(savepath, movfile));
% 
% movfile = sprintf('ten_rpca.%s'     , movformat); im = best_foreground('rpca'     ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results. rpca.    A{best_l_idx('rpca'     )                         }), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_brpca.%s'    , movformat); im = best_foreground('brpca'    ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca'    )}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_irpca_sub.%s', movformat); im = best_foreground('irpca_sub'); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.irpca.sub.A{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_irpca_lin.%s', movformat); im = best_foreground('irpca_lin'); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.irpca.lin.A{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_orpca.%s'    , movformat); im = best_foreground('orpca'    ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca'    )}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_rcpd_sub.%s' , movformat); im = best_foreground('rcpd_sub' ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results. rcpd.sub.A{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_rcpd_lin.%s' , movformat); im = best_foreground('rcpd_lin' ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results. rcpd.lin.A{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}), im), fullfile(savepath, movfile));

%% Latex by Mehdi Bahri
fprintf('\\textbf{Algorithm} & \\textbf{AUC} & \\textbf{Parameter 1} & \\textbf{Parameter 2} & \\textbf{Iterations (indicative)} & \\textbf{Time (indicative)}\\\\ \\hline\n');

fprintf('RPCA2D L1 (Bahri et. al.) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
    best_auc('rpca2d_l1'), params_2d.lambda(best_l_idx('rpca2d_l1')), ...
    results.rpca2d_l1.info{best_l_idx('rpca2d_l1'), best_r_idx('rpca2d_l1')}.iter(end), ...
    results.rpca2d_l1.time(best_l_idx('rpca2d_l1'), best_r_idx('rpca2d_l1')) );

fprintf('Tensor RPCA (CPVR2016) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
    best_auc('tnn'), params_tnn.lambda(best_l_idx('tnn')), ...
    results.tnn.info{best_l_idx('tnn'), best_r_idx('tnn')}.iter(end), ...
    results.tnn.time(best_l_idx('tnn'), best_r_idx('tnn')) );

fprintf('Tensor RPCA (CPVR2014) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
    best_auc('tsvd'), params_tsvd.lambda(best_l_idx('tsvd')), ...
    results.tsvd.info{best_l_idx('tsvd'), best_r_idx('tsvd')}.iter(end), ...
    results.tsvd.time(best_l_idx('tsvd'), best_r_idx('tsvd')) );

fprintf('Tensor HORPCA-S (Goldfarb & Qin) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
    best_auc('horpca_s'), params_horpca.lambda(best_l_idx('horpca_s')), ...
    results.horpca_s.info{best_l_idx('horpca_s'), best_r_idx('horpca_s')}.iter(end), ...
    results.horpca_s.time(best_l_idx('horpca_s'), best_r_idx('horpca_s')) );

fprintf('Cauchy Soft-Thresholding & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
    best_auc('cauchy_st'), params_mest.lambda(best_l_idx('cauchy_st')), params_mest.p(best_r_idx('cauchy_st')), ...
    results.cauchy_st.info{best_l_idx('cauchy_st'), best_r_idx('cauchy_st')}.iter(end), ...
    results.cauchy_st.time(best_l_idx('cauchy_st'), best_r_idx('cauchy_st')) );

fprintf('Welsh Soft-Thresholding & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
    best_auc('welsh_st'), params_mest.lambda(best_l_idx('welsh_st')), params_mest.p(best_r_idx('welsh_st')), ...
    results.welsh_st.info{best_l_idx('welsh_st'), best_r_idx('welsh_st')}.iter(end), ...
    results.welsh_st.time(best_l_idx('welsh_st'), best_r_idx('welsh_st')) );

fprintf('NCTRPCA & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
    best_auc('nctrpca'), params_nctrpca.lambda(best_l_idx('nctrpca')), params_nctrpca.p(best_r_idx('nctrpca')), ...
    results.nctrpca.info{best_l_idx('nctrpca'), best_r_idx('nctrpca')}.iter(end), ...
    results.nctrpca.time(best_l_idx('nctrpca'), best_r_idx('nctrpca')) );