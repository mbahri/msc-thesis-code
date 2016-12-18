%LOAD_FACADE Loads the facade test image and rescales the dynamic range
%
% Mehdi Bahri - Imperial College London
% July, 2016

O = double(imread('facade.png'));
O = normalize_dynamic_range(O);