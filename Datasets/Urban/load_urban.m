%LOAD_URBAN Loads the urban test hyperspectral image and rescales the dynamic range
%
% Mehdi Bahri - Imperial College London
% July, 2016

load('urban.mat');
O = normalize_dynamic_range(urban);
clear urban