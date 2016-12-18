function [X] = add_noise(X, noise_type, param, imidx)
% [X] = add_noise(X, noise_type, param, imidx)
% Adds noise to an array that contains images.
% INPUTS
%       X           the array to add noise to
%       noise_type  the type of noise to add
%       param       noise parameter, specific to the type of noise
%       imidx       the 2 indices that correspond to the images (if needed)
% OUTPUTS
%       X           same as X but with noise added
%
% Georgios Papamakarios
% Imperial College London
% May 2014

switch noise_type
    
    case 'salt_pepper'

        noise_perc = param;
        noise = sign(randn(size(X)));
        noise = (noise + 1) / 2;
        idx = rand(size(X)) > 1 - noise_perc;
        X(idx) = noise(idx);

    case 'patch'
        
        sizeX = size(X);
        imsize = sizeX(imidx);
        X = arr2mat(X, imidx);
        
        for i = 1:size(X,2)
            patch_size = randi(param, [1,2]);
            patch = (sign(randn(patch_size)) + 1) / 2;
            posx = randi(imsize(1) - patch_size(1), 1);
            posy = randi(imsize(2) - patch_size(2), 1);
            im = reshape(X(:,i), imsize);
            im(posx+1:posx+patch_size(1), posy+1:posy+patch_size(2)) = patch;
            X(:,i) = im(:);
        end
        
        X = mat2arr(X, sizeX, imidx);
        
    case 'baboon'
        
        sizeX = size(X);
        imsize = sizeX(imidx);
        X = arr2mat(X, imidx);
        
        patch_size = param;
        patch = double(rgb2gray(imread('baboon.jpeg')));
        patch = patch - min(patch(:));
        patch = patch / max(patch(:));
        patch = imresize(patch, patch_size);
        for i = 1:size(X,2)
            posx = randi(imsize(1) - patch_size(1), 1);
            posy = randi(imsize(2) - patch_size(2), 1);
            im = reshape(X(:,i), imsize);
            im(posx+1:posx+patch_size(1), posy+1:posy+patch_size(2)) = patch;
            X(:,i) = im(:);
        end
        
        X = mat2arr(X, sizeX, imidx);
    
    case 'whole_pics'
        
        sizeX = size(X);
        X = arr2mat(X, imidx);
                
        corrupted_perc = param;
        for i = 1:size(X,2)
            if rand < corrupted_perc 
                X(:,i) = (sign(randn(size(X,1), 1)) + 1) / 2;
            end
        end
        
        X = mat2arr(X, sizeX, imidx);
    
    otherwise
        
        error('Unknown noise type.');
        
end
