clear;
load('/vol/bitbucket/mb2215/experiments/hall/matrices/144x176x300_no_noise.mat')
close all;
% clc;

% choose data
dims = [144 176 300];

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
load_hall
rowidx = [1 2];
colidx = 3;
X = O(:,:,GT_frames);

% Remove unnecessary frames from the results
for i=1:size(results_rpca.A, 1)
    results_rpca.A{i} = results_rpca.A{i}(:,:,GT_frames);
end

imsize = dims(rowidx);
range = [min(X(:)) max(X(:))];

%% best results for each method
best_l_idx = containers.Map;
best_r_idx = containers.Map;
best_tpr = containers.Map;
best_fpr = containers.Map;
best_auc = containers.Map;
best_foreground = containers.Map;

true_mask = GT > min(GT(:)) + 0.1;
threshold = linspace(0, 1, 101);

[best_foreground('rpca'), best_tpr('rpca'), best_fpr('rpca'), best_auc('rpca'), best_l_idx('rpca') ] = select_best_foreground(results_rpca. A, X, true_mask, threshold);