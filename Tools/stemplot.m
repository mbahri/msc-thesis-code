function [ PP ] = stemplot(M, name, lambda)
%STEMPLOT Stem plot the singular values of the input
% saves the plot to a file with name 'fname'

[~, S, ~] = svd(M);
S = diag(S);

PP = stem(S);
title(sprintf('%s | lambda = %g', name, lambda));

% fname = sprintf('STEM_%s_mu_%g_l_%g.pdf', inputname(1), lambda);

% print(fname);
end

