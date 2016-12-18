function [best_X, best_corr, best_i, best_j] = select_best_missing(data, true_data, known_idx, rowidx, colidx)
% [best_X, best_corr, best_i, best_j] = select_best_missing(data, true_data, known_idx, rowidx, colidx)
% Selects the best reconstruction among a 2D array of reconstructions based on the correlation coefficient.
% INPUTS
%       data          2D cell array of reconstructions
%       true_data     original images
%       known_idx     mask of observed pixels
%       rowidx        modes corresponding to the image
%       colidx        modes not corresponding to the image
% OUTPUTS
%       best_X        reconstruction with the highest correlation coefficient
%       best_corr     average correlation coefficient of the best reconstruction
%       best_i        row index of the best foreground
%       best_j        column index of the best foreground
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

best_corr = 0;

A = preprocess_missing(true_data, known_idx, rowidx, colidx);

for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        
        X = preprocess_missing(data{i,j}, known_idx, rowidx, colidx);
        corr = corr_coeff(A, X);
        
        if corr > best_corr
            best_X = X;
            best_corr = corr;
            best_i = i;
            best_j = j;
        end
    end
end 


function [X] = preprocess_missing(data, known_idx, rowidx, colidx)
% Preprocesses an array of images some of which are unobserved.
% INPUTS
%       data          array of images 
%       known_idx     mask of observed pixels
%       rowidx        modes corresponding to the image
%       colidx        modes not corresponding to the image
% OUTPUTS
%       X             a matrix with the preprocessed unobserved images as columns

known_idx = arr2mat(known_idx, rowidx, colidx);
missing_cols = ~logical(sum(known_idx));

X = arr2mat(data, rowidx, colidx);
X = X(:, missing_cols);

for i = 1:size(X,2)
    if median(X(:,i)) < 0, 
        X(:,i) = - X(:,i);
    end
end


function [c] = corr_coeff(X, Y)
% Calculates the average correlation coefficient between the columns of X and Y.
% INPUTS
%       X, Y     matrices of the same size
% OUTPUTS
%       c        average correlation coefficient

N = size(X, 2);

if N ~= size(Y, 2)
    error('Matrices must be of the same size.');
end

all_c = size(1, N);

for i = 1:N
    x = X(:,i);
    y = Y(:,i);
    
    x = x - mean(x);
    y = y - mean(y);
    
    x = x / norm(x);
    y = y / norm(y);
    
    all_c(i) = x' * y;
end

c = mean(all_c);
c = abs(c);
