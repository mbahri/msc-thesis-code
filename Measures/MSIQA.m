function [ out ] = MSIQA( Original, Image )
%MSIQA Wrapper for MSIQA that rescales images to the right dynamic range
% also adds the reconstruction error in Frobenius norm
%
% Mehdi Bahri - Imperial College London
% July, 2016

O1 = rescale(Original, 0, 255);
I1 = rescale(Image, 0, 255);

out = MSIQA_255(O1, I1);

if size(Original, 3) > 1
    out.rel_norm = tensor_relative_error(Original, Image);
else
    out.rel_norm = matrix_relative_error(Original, Image);
end

end

