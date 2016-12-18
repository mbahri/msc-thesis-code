function [tpr, fpr, auc] = roc_analysis(mask, true_mask, threshold)
% [tpr, fpr, auc] = roc_analysis(mask, true_mask, threshold)
% ROC analysis for a list of threshold values.
% INPUTS
%       mask        input mask (greyscale array)
%       true_mask   true mask (binary array)
%       threshold   a list of threshold values
% OUTPUTS
%       tpr         true positive rate [true positives / true]
%       fpr         false positive rate [false positives / false]
%       auc         area under the curve
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

if ~all(size(mask) == size(true_mask))
    error('Masks must be of the same size.');
end

mask = abs(mask);
mask = mask(:);
true_mask = true_mask(:);

N = length(threshold);
tpr = zeros(1, N);
fpr = zeros(1, N);

% calculate tpr and fpr
for i = 1:N
    posv_mask = mask > threshold(i);
    fals_mask = ~true_mask;
    
    tpr(i) = sum(posv_mask & true_mask) / sum(true_mask);
    fpr(i) = sum(posv_mask & fals_mask) / sum(fals_mask);
end

% calculate auc (using a trapezoid approximation)
[sfpr, idx] = sort(fpr);
stpr = tpr(idx);

auc = stpr(1) * sfpr(1) / 2;
for i = 2:N
    auc = auc + (stpr(i-1) + stpr(i)) * (sfpr(i) - sfpr(i-1)) / 2;
end
auc = auc + (stpr(end) + 1) * (1 - sfpr(end)) / 2;
