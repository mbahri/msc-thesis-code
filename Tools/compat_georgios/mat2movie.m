function [movie] = mat2movie(X, imsize, file)
% [movie] = mat2movie(X, imsize, file)
% Converts matrix X into a movie and optionally saves it as AVI.
% INPUTS
%       X         matrix containing vectorised images as columns
%       imsize    the image size
%       file      file name to save movie (optional)
% OUTPUTS
%       movie     a MATLAB movie
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

if nargin < 3
    file = [];
end

N = size(X, 2);
movarr = zeros([imsize 3 N]);

for i = 1:N
    im = reshape(X(:,i), imsize);
    movarr(:,:,1,i) = im;
    movarr(:,:,2,i) = im;
    movarr(:,:,3,i) = im;
end

movie = immovie(movarr);

if ~isempty(file)
    movie2avi(movie, file);
end
    