% test_ten_rcpd.m
% Tests Robust CP Decomposition on tensors.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
close all;

datapath = '../data/';
datafile = 'yaleb10.mat';
%datafile = 'highway.mat';
%datafile = 'multipie10.mat';

load([datapath datafile]);
imsize = size(X);
imsize = imsize(rowidx);
range = [min(X(:)), max(X(:))];

options.verbose = true;
options.method = 'sub';
options.numComp = 30;
%options.p = 1;
%options.q = 1;
options.lambda = 10;

known_mask = rand(size(X)) < 0.2;
options.error_weights = known_mask;
X(~known_mask) = 0;

tic;
Y = ten_rcpd(X, options);
time = toc;

fprintf('Total time = %g min\n', time/60);
X = arr2mat(X,   rowidx, colidx);
A = arr2mat(Y.A, rowidx, colidx);
E = arr2mat(Y.E, rowidx, colidx);
disp_imdata(cat(3, X, A, E), imsize, {'Original', 'Low-rank', 'Sparse'}, [1 3], range);
