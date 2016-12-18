function [ O, X ] = facade_sp( scale, level )
%FACADE_SP Loads the facade image, rescales it, and adds S&P noise
%
% Mehdi Bahri - Imperial College London
% July, 2016

load_facade;

if ~(scale == 1)
    O = imresize(O, scale);
end

rng('default');
s = rng;
rng(12031992);
X = imnoise(O, 'salt & pepper', level);
rng(s);

end

