function [ Y ] = ssim_downsample( X )
%SSIM_DOWNSAMPLE Mean filtering and downsampling for SSIM
% Taken from the author's website:
% Suggested Usage
% 
% The above (ssim_index.m) is a single scale version of the SSIM indexing 
% measure, which is most effective if used at the appropriate scale. 
% The precisely “right” scale depends on both the image resolution and the 
% viewing distance and is usually difficult to be obtained.  In practice, 
% we suggest to use the following empirical formula to determine the scale 
% for images viewed from a typical distance (say 3~5 times of the image 
% height or width): 1) Let F = max(1, round(N/256)), where N is the number 
% of pixels in image height (or width); 2) Average local F by F pixels and 
% then downsample the image by a factor of F; and 3) apply the ssim_index.m
% program. For example, for an 512 by 512 image, F = max(1, round(512/256))
% = 2, so the image should be averaged within a 2 by 2 window and 
% downsampled by a factor of 2 before applying ssim_index.m.
% 
% Mehdi Bahri
% Imperial College London
% August, 2016

N = max(size(X,1), size(X,2));
F = max(1, round(N/256));

h = fspecial('average', [F F]);

Y = imresize(filter2(h, X), 1/F);
end

