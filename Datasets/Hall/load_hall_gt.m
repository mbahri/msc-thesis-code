%LOAD_HALL Loads the hall test video and rescales the dynamic range
%
% Mehdi Bahri - Imperial College London
% July, 2016

load('airport_hall.mat');
O = normalize_dynamic_range(X);
clear X