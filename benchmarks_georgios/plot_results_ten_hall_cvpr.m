% Plots the results of the background subtraction experiment on the hall
% video using tensor methods.
%
% Mehdi Bahri
% Imperial College London
% Aug 2016

clear;
close all;
clc;

% choose data
dims = [144 176 300];

%% Load tensor results

method_type = 'tensors';

% load data
path = '/vol/bitbucket/mb2215/experiments/hall';
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

%% Load baseline

method_type = 'matrices';

% load data
path = '/vol/bitbucket/mb2215/experiments/hall';
file1 = sprintf('%dx%dx%d_no_noise_rnndl.mat', dims);
file2 = sprintf('%dx%dx%d_no_noise_rpca.mat', dims);

load(fullfile(path, method_type, file1));
load(fullfile(path, method_type, file2));

load(fullfile(path, 'params_no_noise_rnndl.mat'));
load(fullfile(path, 'params_no_noise_rpca.mat'));

%% Cleanup stuff from Georgios and others non-present in the paper
results = rmfield(results, 'rpca');
results = rmfield(results, 'brpca');
results = rmfield(results, 'orpca');
results = rmfield(results, 'irpca');
results = rmfield(results, 'rcpd');
results = rmfield(results, 'brtf');
results = rmfield(results, 'rpca2d_l2');
results = rmfield(results, 'rpca2d_gl_l1');
results = rmfield(results, 'rpca2d_gl_l2');

results.rpca = results_rpca;
results.rnndl = results_rnndl;

%%


% plot parameters
axis_fontsize = 14;
linewidth = 2;
markersize_big = 5;
markersize_sml = 3.5;
% colour = distinguishable_colors(9);
colour = parula(9);
% colour = colour([1 5 4 3 7 6 2 8], :);
colour = colour(1:9, :);
marker = {'d', 'x', 's', '>', '<', '+', '*', 'p'};
linestyle = ':';
position = [0 0 23 13];
units = 'centimeters';
paperposition = [0 0 23.5 13.5];

% display pics params
% load('../data/hall.mat');
load_hall
X = O(:,:,GT_frames);
rowidx = [1 2];
colidx = 3;

%% Pre-process the recent matrix baselines
for n={'rpca', 'rnndl'}
    for i=1:size(results.(n{:}).A, 1)
        results.(n{:}).A{i} = results.(n{:}).A{i}(:,:,GT_frames);
    end
end

imsize = dims(rowidx);
range = [min(X(:)) max(X(:))];

% figure save path
savepath = fullfile('/vol/bitbucket/mb2215/Thesis', 'results/cvpr', 'hall');

%% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_tpr   = containers.Map;
best_fpr   = containers.Map;
best_auc   = containers.Map;
best_foreground = containers.Map;

true_mask = GT > min(GT(:)) + 0.1;
threshold = linspace(0, 1, 101);

% Matrix baselines
[best_foreground('rpca'     ), best_tpr('rpca'     ), best_fpr('rpca'     ), best_auc('rpca'     ), best_l_idx('rpca'     ), best_r_idx('rpca'    )] = select_best_foreground(results. rpca.    A, X, true_mask, threshold);
[best_foreground('rnndl'    ), best_tpr('rnndl'    ), best_fpr('rnndl'    ), best_auc('rnndl'    ), best_l_idx('rnndl'    ), best_r_idx('rnndl'    )] = select_best_foreground(results.rnndl.    A, X, true_mask, threshold);

% Me
[best_foreground('rpca2d_l1' ), best_tpr('rpca2d_l1' ), best_fpr('rpca2d_l1' ), best_auc('rpca2d_l1' ), best_l_idx('rpca2d_l1' ), best_r_idx('rpca2d_l1' )] = select_best_foreground(results.rpca2d_l1.A, X, true_mask, threshold);

% The rest
[best_foreground('tnn' ), best_tpr('tnn' ), best_fpr('tnn' ), best_auc('tnn' ), best_l_idx('tnn' ), best_r_idx('tnn' )] = select_best_foreground(results.tnn.A, X, true_mask, threshold);
[best_foreground('tsvd' ), best_tpr('tsvd' ), best_fpr('tsvd' ), best_auc('tsvd' ), best_l_idx('tsvd' ), best_r_idx('tsvd' )] = select_best_foreground(results.tsvd.A, X, true_mask, threshold);

[best_foreground('cauchy_st' ), best_tpr('cauchy_st' ), best_fpr('cauchy_st' ), best_auc('cauchy_st' ), best_l_idx('cauchy_st' ), best_r_idx('cauchy_st' )] = select_best_foreground(results.cauchy_st.A, X, true_mask, threshold);
[best_foreground('welsh_st' ), best_tpr('welsh_st' ), best_fpr('welsh_st' ), best_auc('welsh_st' ), best_l_idx('welsh_st' ), best_r_idx('welsh_st' )] = select_best_foreground(results.welsh_st.A, X, true_mask, threshold);
[best_foreground('horpca_s' ), best_tpr('horpca_s' ), best_fpr('horpca_s' ), best_auc('horpca_s' ), best_l_idx('horpca_s' ), best_r_idx('horpca_s' )] = select_best_foreground(results.horpca_s.A, X, true_mask, threshold);

[best_foreground('nctrpca' ), best_tpr('nctrpca' ), best_fpr('nctrpca' ), best_auc('nctrpca' ), best_l_idx('nctrpca' ), best_r_idx('nctrpca' )] = select_best_foreground(results.nctrpca.A, X, true_mask, threshold);

maxiter = max([
               results. rpca2d_l1.info{best_l_idx('rpca2d_l1' ), best_r_idx('rpca2d_l1' )}.iter(end) ...
               results. tnn.info{best_l_idx('tnn' ), best_r_idx('tnn' )}.iter(end) ...
               results. horpca_s.info{best_l_idx('horpca_s' ), best_r_idx('horpca_s' )}.iter(end) ...
               results. nctrpca.info{best_l_idx('nctrpca' ), best_r_idx('nctrpca' )}.iter(end) ...
               results. rpca.info{best_l_idx('rpca' ), best_r_idx('rpca' )}.iter(end) ...
               results. rnndl.info{best_l_idx('rnndl' ), best_r_idx('rnndl' )}.iter(end) ...
               ] );


save('roc_analysis_hall',  ...
    'best_l_idx', ...
    'best_r_idx', ...
    'best_tpr', ...
    'best_fpr', ...
    'best_auc', ...
    'best_foreground', ...
    '-v7.3' ...
);

%% ROC curves
cnt = 1;

fig = figure; 
fpr = norminv(best_fpr('rpca'     )); fnr = norminv(1 - best_tpr('rpca'     )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; hold on;
fpr = norminv(best_fpr('rnndl'    )); fnr = norminv(1 - best_tpr('rnndl'    )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{1}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('rpca2d_l1' )); fnr = norminv(1 - best_tpr('rpca2d_l1' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{2}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1;  hold on;

fpr = norminv(best_fpr('tnn' )); fnr = norminv(1 - best_tpr('tnn' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{3}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('tsvd' )); fnr = norminv(1 - best_tpr('tsvd' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{4}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('horpca_s' )); fnr = norminv(1 - best_tpr('horpca_s' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{5}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('cauchy_st' )); fnr = norminv(1 - best_tpr('cauchy_st' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{7}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 
fpr = norminv(best_fpr('welsh_st' )); fnr = norminv(1 - best_tpr('welsh_st' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{7}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth); cnt = cnt + 1; 

fpr = norminv(best_fpr('nctrpca' )); fnr = norminv(1 - best_tpr('nctrpca' )); idx = (abs(fpr) ~= inf) & (abs(fnr) ~= inf); plot(fpr(idx), fnr(idx), linestyle, 'Marker', marker{8}, 'Color', colour(cnt,:), 'MarkerSize', markersize_sml, 'LineWidth', linewidth);

ylabel('False Negative Rate', 'FontSize', axis_fontsize);
xlabel('False Positive Rate', 'FontSize', axis_fontsize);
legend('RPCA (baseline)', 'RNNDL (baseline)', 'KDRSDL (proposed)', 'CVPR 2016', 'CVPR 2014', ...
    'HORPCA-S', 'Cauchy ST', 'Welsh ST', 'NC TRPCA', ...
    'Location', 'Best');

set(gca, 'FontSize', axis_fontsize);
ppxtick = @(x) ( strrep(sprintf('%.0e',x), 'e-0', 'e-') );
xtick = [10.^(-4:-1) 0.5]; set(gca, 'XTick', norminv(xtick)); set(gca, 'XTickLabel', arrayfun(ppxtick, xtick, 'UniformOutput', 0));
ytick = 0.1 : 0.2 : 0.9;   set(gca, 'YTick', norminv(ytick)); set(gca, 'YTickLabel', {ytick});
xlim(norminv([2.2421e-05 0.9846])); % xlim(norminv([2.2421e-05 0.8318]));
ylim(norminv([0.0012 0.9998])); % ylim(norminv([0.0364 0.9998]));
set(fig, 'PaperPositionMode', 'auto'); 
% saveas(fig, fullfile(savepath, 'ten_hall_roc_curves.eps'), 'epsc');
% saveas(fig, fullfile(savepath, 'ten_hall_roc_curves_nb.eps'), 'eps');
saveas(fig, fullfile(savepath, 'ten_hall_roc_curves_parula.eps'), 'epsc');
saveas(fig, fullfile(savepath, 'ten_hall_roc_curves_parula_nb.eps'), 'eps');

fig.PaperOrientation = 'landscape';
fig.Renderer = 'painters';
fig.Units = units;
fig.Position = position;
fig.PaperUnits = units;
fig.PaperPosition = paperposition;
print(fig, '-depsc', fullfile(savepath, 'ten_hall_roc_curves_parula'));
print(fig, '-deps', fullfile(savepath, 'ten_hall_roc_curves_parula_nb'));
% print(fig, '-depsc', fullfile(savepath, 'ten_hall_roc_curves'));
% print(fig, '-deps', fullfile(savepath, 'ten_hall_roc_curves_nb'));

%% best results

fprintf('RPCA2D L1 (Bahri et. al.): \n');
fprintf('\t lambda = %g \n', params_2d.lambda(best_l_idx('rpca2d_l1')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('rpca2d_l1')));
fprintf('\t time = %g sec \n', results.rpca2d_l1.time(best_l_idx('rpca2d_l1'),best_r_idx('rpca2d_l1')));
fprintf('\t auc = %g \n', best_auc('rpca2d_l1'));
fprintf('\t iterations = %d \n\n', results.rpca2d_l1.info{best_l_idx('rpca2d_l1'),best_r_idx('rpca2d_l1')}.iter(end));

fprintf('Tensor RPCA (CPVR 2016): \n');
fprintf('\t lambda = %g \n', params_tnn.lambda(best_l_idx('tnn')));
% fprintf('\t rank = %g \n', ten_rank(best_r_idx('tnn')));
fprintf('\t time = %g sec \n', results.tnn.time(best_l_idx('tnn'),best_r_idx('tnn')));
fprintf('\t auc = %g \n', best_auc('tnn'));
fprintf('\t iterations = %d \n\n', results.tnn.info{best_l_idx('tnn'),best_r_idx('tnn')}.iter(end));

fprintf('Tensor RPCA (CPVR 2014): \n');
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

fprintf('RPCA: \n');
fprintf('\t threshold = %g \n', params_rpca.p(best_r_idx('rpca')));
fprintf('\t rank = %g \n', params_rpca.lambda(best_l_idx('rpca')));
fprintf('\t time = %g sec \n', results.rpca.time(best_l_idx('rpca'),best_r_idx('rpca')));
fprintf('\t auc = %g \n', best_auc('rpca'));
fprintf('\t iterations = %d \n\n', results.rpca.info{best_l_idx('rpca'),best_r_idx('rpca')}.iter(end));

fprintf('RNNDL: \n');
fprintf('\t threshold = %g \n', params_rnndl.p(best_r_idx('rnndl')));
fprintf('\t rank = %g \n', params_rnndl.lambda(best_l_idx('rnndl')));
fprintf('\t time = %g sec \n', results.rnndl.time(best_l_idx('rnndl'),best_r_idx('rnndl')));
fprintf('\t auc = %g \n', best_auc('rnndl'));
fprintf('\t iterations = %d \n\n', results.rnndl.info{best_l_idx('rnndl'),best_r_idx('rnndl')}.iter(end));

%% save images
algorithms = fieldnames(results);
im_idx = [1 2];
for nn=algorithms'
    n = nn{:}
    for i=im_idx
        imfile = sprintf('hall_bg_%d_%s.png', i, n);
        im = results.(n).A{best_l_idx(n), best_r_idx(n)}(:,:,i);
        imwrite(double(im), fullfile(savepath, imfile), 'png');
        
        imfile = sprintf('hall_fg_%d_%s.png', i, n);
        im = best_foreground(n);
        im = im(:,:,i);
        imwrite(double(im), fullfile(savepath, imfile), 'png');
    end
end

fprintf('End\n');

%% save images
% im_idx = [130 164 315 355];
% imformat = 'png';
% 
% for i = 1:length(im_idx)
%     imfile = sprintf('ten_hall_original_%d.%s', im_idx(i), imformat); im = X        (:,:,im_idx(i)); imwrite(im, fullfile(savepath, imfile), imformat);
%     imfile = sprintf('ten_hall_mask_%d.%s'    , im_idx(i), imformat); im = true_mask(:,:,im_idx(i)); imwrite(im, fullfile(savepath, imfile), imformat);
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
% movfile = sprintf('ten_hall.%s', movformat); arr2movie(cat(2, X, true_mask), fullfile(savepath, movfile));
% 
% movfile = sprintf('ten_rpca.%s'     , movformat); im = best_foreground('rpca'     ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results. rpca.    A{best_l_idx('rpca'     )                         }), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_brpca.%s'    , movformat); im = best_foreground('brpca'    ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.brpca.    A{best_l_idx('brpca'    ), best_r_idx('brpca'    )}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_irpca_sub.%s', movformat); im = best_foreground('irpca_sub'); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.irpca.sub.A{best_l_idx('irpca_sub'), best_r_idx('irpca_sub')}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_irpca_lin.%s', movformat); im = best_foreground('irpca_lin'); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.irpca.lin.A{best_l_idx('irpca_lin'), best_r_idx('irpca_lin')}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_orpca.%s'    , movformat); im = best_foreground('orpca'    ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results.orpca.    A{best_l_idx('orpca'    ), best_r_idx('orpca'    )}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_rcpd_sub.%s' , movformat); im = best_foreground('rcpd_sub' ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results. rcpd.sub.A{best_l_idx('rcpd_sub' ), best_r_idx('rcpd_sub' )}), im), fullfile(savepath, movfile));
% movfile = sprintf('ten_rcpd_lin.%s' , movformat); im = best_foreground('rcpd_lin' ); im = im - min(im(:)); im = im / max(im(:)); arr2movie(cat(2, double(results. rcpd.lin.A{best_l_idx('rcpd_lin' ), best_r_idx('rcpd_lin' )}), im), fullfile(savepath, movfile));

% %% Latex by Mehdi Bahri
% fprintf('\\textbf{Algorithm} & \\textbf{AUC} & \\textbf{Parameter 1} & \\textbf{Parameter 2} & \\textbf{Iterations (indicative)} & \\textbf{Time (indicative)}\\\\ \\hline\n');
% 
% % fprintf('RPCA2D L2 (Bahri et. al.) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
% %     best_auc('rpca2d_l2'), params_2d.lambda(best_l_idx('rpca2d_l2')), ...
% %     results.rpca2d_l2.info{best_l_idx('rpca2d_l2'), best_r_idx('rpca2d_l2')}.iter(end), ...
% %     results.rpca2d_l2.time(best_l_idx('rpca2d_l2'), best_r_idx('rpca2d_l2')) );
% 
% fprintf('RPCA2D L1 (Bahri et. al.) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('rpca2d_l1'), params_2d.lambda(best_l_idx('rpca2d_l1')), ...
%     results.rpca2d_l1.info{best_l_idx('rpca2d_l1'), best_r_idx('rpca2d_l1')}.iter(end), ...
%     results.rpca2d_l1.time(best_l_idx('rpca2d_l1'), best_r_idx('rpca2d_l1')) );
% 
% % fprintf('RPCA2D GL L2 (Bahri et. al.) & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
% %     best_auc('rpca2d_gl_l2'), params_2d.lambda(best_l_idx('rpca2d_gl_l2')), params_2d.p(best_r_idx('rpca2d_gl_l2')), ...
% %     results.rpca2d_gl_l2.info{best_l_idx('rpca2d_gl_l2'), best_r_idx('rpca2d_gl_l2')}.iter(end), ...
% %     results.rpca2d_gl_l2.time(best_l_idx('rpca2d_gl_l2'), best_r_idx('rpca2d_gl_l2')) );
% % 
% % fprintf('RPCA2D GL L1 (Bahri et. al.) & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
% %     best_auc('rpca2d_gl_l1'), params_2d.lambda(best_l_idx('rpca2d_gl_l1')), params_2d.p(best_r_idx('rpca2d_gl_l1')), ...
% %     results.rpca2d_gl_l1.info{best_l_idx('rpca2d_gl_l1'), best_r_idx('rpca2d_gl_l1')}.iter(end), ...
% %     results.rpca2d_gl_l1.time(best_l_idx('rpca2d_gl_l1'), best_r_idx('rpca2d_gl_l1')) );
% 
% fprintf('Tensor RPCA (CPVR2016) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('tnn'), params_tnn.lambda(best_l_idx('tnn')), ...
%     results.tnn.info{best_l_idx('tnn'), best_r_idx('tnn')}.iter(end), ...
%     results.tnn.time(best_l_idx('tnn'), best_r_idx('tnn')) );
% 
% fprintf('Tensor RPCA (CPVR2014) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('tsvd'), params_tsvd.lambda(best_l_idx('tsvd')), ...
%     results.tsvd.info{best_l_idx('tsvd'), best_r_idx('tsvd')}.iter(end), ...
%     results.tsvd.time(best_l_idx('tsvd'), best_r_idx('tsvd')) );
% 
% fprintf('Tensor HORPCA-S (Goldfarb & Qin) & %.4f & %.3g & & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('horpca_s'), params_horpca.lambda(best_l_idx('horpca_s')), ...
%     results.horpca_s.info{best_l_idx('horpca_s'), best_r_idx('horpca_s')}.iter(end), ...
%     results.horpca_s.time(best_l_idx('horpca_s'), best_r_idx('horpca_s')) );
% 
% fprintf('Cauchy Soft-Thresholding & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('cauchy_st'), params_mest.lambda(best_l_idx('cauchy_st')), params_mest.p(best_r_idx('cauchy_st')), ...
%     results.cauchy_st.info{best_l_idx('cauchy_st'), best_r_idx('cauchy_st')}.iter(end), ...
%     results.cauchy_st.time(best_l_idx('cauchy_st'), best_r_idx('cauchy_st')) );
% 
% fprintf('Welsh Soft-Thresholding & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('welsh_st'), params_mest.lambda(best_l_idx('welsh_st')), params_mest.p(best_r_idx('welsh_st')), ...
%     results.welsh_st.info{best_l_idx('welsh_st'), best_r_idx('welsh_st')}.iter(end), ...
%     results.welsh_st.time(best_l_idx('welsh_st'), best_r_idx('welsh_st')) );
% 
% fprintf('NCTRPCA & %.4f & %.3g & %.3g & %u & %.2f\\\\ \\hline\n', ...
%     best_auc('nctrpca'), params_nctrpca.lambda(best_l_idx('nctrpca')), params_nctrpca.p(best_r_idx('nctrpca')), ...
%     results.nctrpca.info{best_l_idx('nctrpca'), best_r_idx('nctrpca')}.iter(end), ...
%     results.nctrpca.time(best_l_idx('nctrpca'), best_r_idx('nctrpca')) );
