function [ output_args ] = show_3( A, B, C )
%SHOW_3 Shows three images side by side
%
% Mehdi Bahri - Imperial College London
% July, 2016

subplot(1,3,1), ddisp(A)
subplot(1,3,2), ddisp(B)
subplot(1,3,3), ddisp(C)

% subplot(1,3,1), imshow(A, [0, 1])
% subplot(1,3,2), imshow(B, [0, 1])
% subplot(1,3,3), imshow(C, [0, 1])

end

