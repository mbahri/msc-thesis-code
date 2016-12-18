%LOAD_HIGHWAY1 Loads the basketball test video and rescales the dynamic range
%
% Mehdi Bahri - Imperial College London
% July, 2016

load('highway_video.mat');
O = normalize_dynamic_range(T);
clear T
