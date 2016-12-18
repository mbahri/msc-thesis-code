function [ O, X ] = yale_patch( scale, max_size )
%YALE_PATCH Loads the Yale database, rescales it, and adds random patchs of
%maximum size max_size
%
% Mehdi Bahri - Imperial College London
% July, 2016

load_yale;

if ~(scale == 1)
    O = imresize(O, scale);
end

rng('default');
s = rng;
rng(12031992);
X = add_noise(O, 'patch', max_size, [1 2]);
rng(s);

end

