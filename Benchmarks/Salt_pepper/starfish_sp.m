function [ O, X ] = starfish_sp( scale, level )
%STARFISH_SP Loads the starfish image, rescales it, and adds S&P noise
%
% Mehdi Bahri - Imperial College London
% July, 2016

load_starfish;

if ~(scale == 1)
    O = imresize(O, scale);
end

rng('default');
s = rng;
rng(12031992);
X = imnoise(O, 'salt & pepper', level);
rng(s);

end

