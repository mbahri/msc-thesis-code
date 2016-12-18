function [movie] = arr2movie(X, file)
% [movie] = arr2movie(X, file)
% Converts array X into a movie and optionally saves it as AVI.
% INPUTS
%       X         array containing images as frontal slices
%       file      file name to save movie (optional)
% OUTPUTS
%       movie     a MATLAB movie
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

if nargin < 2
    file = [];
end

sizX = size(X);
imsize = sizX(1:2);
N = sizX(3);
movarr = zeros([imsize 3 N]);

for i = 1:N
    im = X(:,:,i);
    movarr(:,:,1,i) = im;
    movarr(:,:,2,i) = im;
    movarr(:,:,3,i) = im;
end

movie = immovie(movarr);

if ~isempty(file)
    movie2avi(movie, file);
end
    