function PSNR = PSNR(origImg, distImg, max_o )
%PSNR Computes the Peak Signal to Noise ratio.
% Based on an existing code. Modified by:
%
% Mehdi Bahri - Imperial College London
% July, 2016

if nargin<3
    max_o = 255;
end

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;
MSE = sum(sum(error .* error)) / (M * N);

if(MSE > 0)
    PSNR = 10*log10(max_o*max_o/MSE);
else
    PSNR = 99;
end
