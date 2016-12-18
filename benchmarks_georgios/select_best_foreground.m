function [best_foreground, best_tpr, best_fpr, best_auc, best_i, best_j] = select_best_foreground(background, origdata, true_mask, threshold)
% [best_foreground, best_tpr, best_fpr, best_auc, best_i, best_j] = select_best_foreground(background, origdata, true_mask, threshold)
% Selects the best foreground among a 2D array of backgrounds based on ROC analysis.
% INPUTS
%       background         2D cell array of backgrounds (greyscale arrays)
%       origdata           original array of both foreground and background
%       true_mask          true mask (binary array)
%       threshold          a list of threshold values for the ROC analysis
% OUTPUTS
%       best_foreground    foreground with the largest area under the curve
%       best_tpr           true positive rate of the best foreground
%       best_fpr           false positive rate of the best foreground
%       best_auc           area under the curve of the best foreground
%       best_i             row index of the best foreground
%       best_j             column index of the best foreground
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

best_auc = -inf;

for i = 1:size(background, 1)
    for j = 1:size(background, 2)
        
        foreground = origdata - background{i,j};
        foreground = double(foreground);
        [tpr, fpr, auc] = roc_analysis(foreground, true_mask, threshold);
        
        if auc > best_auc
            best_foreground = foreground;
            best_tpr = tpr;
            best_fpr = fpr;
            best_auc = auc;
            best_i = i;
            best_j = j;
        end
    end
end            
