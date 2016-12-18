function [ output_args ] = show_4( A, B, C, D )
%SHOW_4 Show four images side by side
%
% Mehdi Bahri - Imperial College London
% July, 2016

subplot(1,4,1), ddisp(A);
subplot(1,4,2), ddisp(B);
subplot(1,4,3), ddisp(C);
subplot(1,4,4), ddisp(D);

end

