%LOAD_STARFISH Loads the starfish test image and rescales the dynamic range
%
% Mehdi Bahri - Imperial College London
% July, 2016

O = double(imread('starfish.jpg'));
O = normalize_dynamic_range(O);
