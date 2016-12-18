clear
close all

subfolder = 'Yale_SP';
dataset = 'yale';
noise_type = 'sp';
noise_level = 0.3;

% Load files, extract only the useful information, clear the rest to avoid going OOM

fprintf('Loading BRTF...\n');
method = 'brtf';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures.brtf, best.brtf] = distilate(results_brtf, 'lambda');
clear results_brtf;

fprintf('Loading Cauchy ST...\n');
method = 'cauchy_st';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)

mask = ~cellfun(@isempty, results_cauchy_st(:));

temp = results_cauchy_st(mask);
temp = reshape(temp, 20, 20)';

results_cauchy_st = temp;

[measures_2d.cauchy_st, best_2d.cauchy_st] = distilate_2d(results_cauchy_st, 'sigma', 'alpha');

save(to_load, 'results_cauchy_st', '-v7.3');

clear results_cauchy_st temp mask;

fprintf('Loading Georgios...\n');
method = 'georgios';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures.ten_rpca, best.ten_rpca] = distilate(results_georgios.ten_rpca, 'lambda');
[measures.ten_brpca, best.ten_brpca] = distilate(results_georgios.ten_brpca, 'lambda');
[measures.ten_orpca, best.ten_orpca] = distilate(results_georgios.ten_orpca, 'lambda');
[measures.ten_rcpd, best.ten_rcpd] = distilate(results_georgios.ten_rcpd, 'lambda');
clear results_georgios;

fprintf('Loading HORPCA-S...\n');
method = 'horpca_s';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures.horpca_s, best.horpca_s] = distilate(results_horpca_s.horpca_s, 'lambda');
[measures.horpca_s_tc, best.horpca_s_tc] = distilate(results_horpca_s.horpca_s_tc, 'lambda');
clear results_horpca_s

fprintf('Loading Mehdi...\n');
method = 'mehdi';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures.rpca2d_l1, best.rpca2d_l1] = distilate(results_mehdi.rpca2d_l1, 'lambda');
[measures.rpca2d_l2, best.rpca2d_l2] = distilate(results_mehdi.rpca2d_l2, 'lambda');
[measures.rpca2d_gl_l1, best.rpca2d_gl_l1] = distilate(results_mehdi.rpca2d_gl_l1, 'lambda');
[measures.rpca2d_gl_l2, best.rpca2d_gl_l2] = distilate(results_mehdi.rpca2d_gl_l2, 'lambda');
clear results_mehdi;

fprintf('Loading CVPR2016...\n');
method = 'tnn';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures.cvpr2016_tnn, best.cvpr2016_tnn] = distilate(results_trpca_tnn, 'lambda');
clear results_trpca_tnn;

fprintf('Loading CVPR2014...\n');
method = 'tsvd';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures.cvpr2014_tsvd, best.cvpr2014_tsvd] = distilate(results_trpca_tsvd, 'lambda');
clear results_trpca_tsvd;

% measures.cvpr2014_tsvd.param1 = zeros(1,40);
% measures.cvpr2014_tsvd.psnr = zeros(1,40);
% measures.cvpr2014_tsvd.fsim = zeros(1,40);
% measures.cvpr2014_tsvd.ssim = zeros(1,40);
% measures.cvpr2014_tsvd.msam = ones(1,40);
% measures.cvpr2014_tsvd.rel_norm = ones(1,40);
% 
% best.cvpr2014_tsvd.psnr = struct('value', 0, 'index', 1, 'param1',  2.2204e-16);
% best.cvpr2014_tsvd.fsim = struct('value', 0, 'index', 1, 'param1',  2.2204e-16);
% best.cvpr2014_tsvd.ssim = struct('value', 0, 'index', 1, 'param1',  2.2204e-16);
% best.cvpr2014_tsvd.msam = struct('value', 0, 'index', 1, 'param1',  1);
% best.cvpr2014_tsvd.rel_norm = struct('value', 0, 'index', 1, 'param1',  1);


fprintf('Loading Welsh ST...\n');
method = 'welsh_st';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures_2d.welsh_st, best_2d.welsh_st] = distilate_2d(results_welsh_st, 'sigma', 'alpha');
clear results_welsh_st;

fprintf('NC TRPCA...\n');
method = 'nctrpca';
to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
load(to_load)
[measures_2d.nctrpca, best_2d.nctrpca] = distilate_2d(results, 'threshold', 'rank');
clear results;

% measures_2d = struct();
% best_2d = struct();

save(sprintf('preprocess_%s_%s_%f.mat', dataset, noise_type, noise_level), 'measures', 'measures_2d', 'best', 'best_2d', '-v7.3');

%load 'loaded_results_facade_preprocess.mat'

algorithms = {'ten_rpca',  'ten_brpca', 'ten_orpca', 'ten_rcpd', ...
              'rpca2d_l1', 'rpca2d_l2', 'rpca2d_gl_l1', 'rpca2d_gl_l2', ...
              'cvpr2016_tnn', 'brtf', 'horpca_s', 'horpca_s_tc', 'cvpr2014_tsvd'};
algorithms_2d = {'welsh_st', 'cauchy_st', 'nctrpca'};

n_1d = length(algorithms);
n_2d = length(algorithms_2d);

legend_names = {'Tensor RPCA (Papamakarios et. al.)', ...
                'Tensor BRPCA (Papamakarios et. al.)', ...
                'Tensor RHOSVD (Papamakarios et. al.)', ...
                'Tensor RCPD (Papamakarios et. al.)', ...
                '2D RPCA L1 (Bahri et. al.)', ...
                '2D RPCA L2 (Bahri et. al.)', ...
                '2D RPCA GL L1 (Bahri et. al.)', ...
                '2D RPCA GL L2 (Bahri et. al.)', ...
                'Tensor RPCA (CVPR 2016)', ...
                'BRTF', ...
                'HORPCA-S', ...
                'HORPCA-S TC', ...
                'Tensor RPCA (CVPR 2014)', ...
                'Welsh ST', ...
                'Cauchy ST', ...
                'NC TRPCA' ...
                };
legend_symbols = {'*', '*', '*', '*', '+', '+', '+', '+', 'o', 'x', 'v', 'v', 's', 'd', 'd', '^'};
legend_colours = distinguishable_colors(length(legend_names));

% subplot(1,3,1);
% figure
% plot(lambda, ten_brpca_psnr, '-o', lambda, ten_rpca_psnr, '-o', lambda, ten_rcpd_psnr, '-o', lambda, ten_orpca_psnr, '-o', ...
%     lambda, rpca2d_gl_l1_psnr, '-+', lambda, rpca2d_l2_psnr, '-+', lambda, rpca2d_l1_psnr, '-+', lambda, rpca2d_gl_l2_psnr, '-+', ...
%     lambda_tnn, trpca_tnn_psnr, '-x');
% legend('brpca', 'rpca', 'rcpd', 'orpca', '2D GL L1', '2D L2', '2D L1', '2D GL L2', 'TRPCA TNN');
% title('PSNR (more is better)');
% xlim([0 6e-2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSNR
fig_psnr = figure
for i=1:n_1d
    n = algorithms{i};
    semilogx(measures.(n).param1, measures.(n).psnr, 'Marker', legend_symbols{i}, 'Color', legend_colours(i, :)); hold on;
end
for i=1:n_2d
    j = n_1d + i;
    n = algorithms_2d{i};
    semilogx(measures_2d.(n).param2, measures_2d.(n).psnr(best_2d.(n).psnr.index_param1, :), 'Marker', legend_symbols{j}, 'Color', legend_colours(j, :)); hold on;
end
title('PSNR');
xlim([2e-4, 5e-1]);
legend(legend_names{:}, 'location', 'northeastoutside');
fig_psnr.Units = 'centimeters';
fig_psnr.Position = [0 0 29.7 21];
fig_psnr.PaperPositionMode = 'auto';
% set(fig_psnr, 'PaperPositionMode', 'auto'); 
fig_psnr.PaperOrientation = 'landscape';
fig_psnr.Renderer = 'painters';
fig_psnr.PaperUnits = 'centimeters';
fig_psnr.PaperPosition = [0 0 29.7 21];
file = sprintf('psnr_%s_%s_%f.pdf', dataset, noise_type, noise_level);
print(fig_psnr, '-dpdf', file);
% saveas(fig_psnr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSIM
fig_ssim = figure
for i=1:n_1d
    n = algorithms{i};
    semilogx(measures.(n).param1, measures.(n).ssim, 'Marker', legend_symbols{i}, 'Color', legend_colours(i, :)); hold on;
end
for i=1:n_2d
    j = n_1d + i;
    n = algorithms_2d{i};
    semilogx(measures_2d.(n).param2, measures_2d.(n).ssim(best_2d.(n).ssim.index_param1, :), 'Marker', legend_symbols{j}, 'Color', legend_colours(j, :)); hold on;
end
title('SSIM');
xlim([2e-4, 5e-1]);
legend(legend_names{:}, 'location', 'northeastoutside');
fig_ssim.Units = 'centimeters';
fig_ssim.Position = [0 0 29.7 21];
fig_ssim.PaperPositionMode = 'auto';
fig_ssim.PaperOrientation = 'landscape';
fig_ssim.Renderer = 'painters';
fig_ssim.PaperUnits = 'centimeters';
fig_ssim.PaperPosition = [0 0 29.7 21];
file = sprintf('ssim_%s_%s_%f.pdf', dataset, noise_type, noise_level);
print(fig_ssim, '-dpdf', file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FSIM
fig_fsim = figure
for i=1:n_1d
    n = algorithms{i};
    semilogx(measures.(n).param1, measures.(n).fsim, 'Marker', legend_symbols{i}, 'Color', legend_colours(i, :)); hold on;
end
for i=1:n_2d
    j = n_1d + i;
    n = algorithms_2d{i};
    semilogx(measures_2d.(n).param2, measures_2d.(n).fsim(best_2d.(n).fsim.index_param1, :), 'Marker', legend_symbols{j}, 'Color', legend_colours(j, :)); hold on;
end
title('FSIM');
xlim([2e-4, 5e-1]);
legend(legend_names{:}, 'location', 'northeastoutside');
fig_fsim.Units = 'centimeters';
fig_fsim.Position = [0 0 29.7 21];
fig_fsim.PaperPositionMode = 'auto';
fig_fsim.PaperOrientation = 'landscape';
fig_fsim.Renderer = 'painters';
fig_fsim.PaperUnits = 'centimeters';
fig_fsim.PaperPosition = [0 0 29.7 21];
file = sprintf('fsim_%s_%s_%f.pdf', dataset, noise_type, noise_level);
print(fig_fsim, '-dpdf', file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSAM
fig_msam = figure
for i=1:n_1d
    n = algorithms{i};
    semilogx(measures.(n).param1, measures.(n).msam, 'Marker', legend_symbols{i}, 'Color', legend_colours(i, :)); hold on;
end
for i=1:n_2d
    j = n_1d + i;
    n = algorithms_2d{i};
    semilogx(measures_2d.(n).param2, measures_2d.(n).msam(best_2d.(n).msam.index_param1, :), 'Marker', legend_symbols{j}, 'Color', legend_colours(j, :)); hold on;
end
title('MSAM');
xlim([2e-4, 5e-1]);
ylim([0 6e-1]);
legend(legend_names{:}, 'location', 'northeastoutside');
fig_msam.Units = 'centimeters';
fig_msam.Position = [0 0 29.7 21];
fig_msam.PaperPositionMode = 'auto';
fig_msam.PaperOrientation = 'landscape';
fig_msam.Renderer = 'painters';
fig_msam.PaperUnits = 'centimeters';
fig_msam.PaperPosition = [0 0 29.7 21];
file = sprintf('msam_%s_%s_%f.pdf', dataset, noise_type, noise_level);
print(fig_msam, '-dpdf', file);

fprintf('"****** BEST RESULTS %s %s %f ******"\n', dataset, noise_type, noise_level);

for i=1:n_1d
    n = algorithms{i};
    fprintf('"%s"\tPSNR: %f for param1=%f param2=NA\tSSIM: %f for param1=%f param2=NA\tFSIM: %f for param1=%f param2=NA\tMSAM: %f for param1=%f param2=NA\n', ...
    legend_names{i}, best.(n).psnr.value, best.(n).psnr.param1, best.(n).ssim.value, best.(n).ssim.param1, ...
    best.(n).fsim.value, best.(n).fsim.param1, best.(n).msam.value, best.(n).msam.param1);
end
for i=1:n_2d
    j = n_1d + i;
    n = algorithms_2d{i};
    fprintf('"%s"\tPSNR: %f for param1=%f param2=%f\tSSIM: %f for param1=%f param2=%f\tFSIM: %f for param1=%f param2=%f\tMSAM: %f for param1=%f param2=%f\n', ...
    legend_names{j}, best_2d.(n).psnr.value, best_2d.(n).psnr.param1, best_2d.(n).psnr.param2, ...
    best_2d.(n).ssim.value, best_2d.(n).ssim.param1, best_2d.(n).ssim.param2, ...
    best_2d.(n).fsim.value, best_2d.(n).fsim.param1, best_2d.(n).fsim.param2, ...
    best_2d.(n).msam.value, best_2d.(n).msam.param1, best_2d.(n).msam.param2);
end
