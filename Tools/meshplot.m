function [ PP ] = meshplot(M, mu, lambda)
%Mesh plot the input
% saves the plot to a file with name 'fname'

PP = mesh(M);
title(sprintf('%s | mu = %g | l = %g', inputname(1), mu, lambda));
view(7, 14);

fname = sprintf('MESH_%s_mu_%g_l_%g.pdf', inputname(1), mu, lambda);

% print(fname);
end

