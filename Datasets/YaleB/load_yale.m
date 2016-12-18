%LOAD_YALE Loads the Yale database and cleans unnecessary variables
%
% Mehdi Bahri - Imperial College London
% July, 2016

load('yaleb10_full_res.mat');
O = X;
clear X cids colidx rowidx