function [ C ] = left_right_matrix_tensor_prod( A, B, Tensor, mod_l, mod_r )
%LEFT_RIGHT_MATRIX_TENSOR_PROD L/R product of tensor slices and matrices
%   Multiplies each slice of Tensor by A on the left and B on the right
%   mod_l: transposition modififers for the left product
%   mod_r: transposition modifiers for the right product
%   REQUIRES MMX

C = mmx('mult', A, Tensor, mod_l);
C = mmx('mult', C, B, mod_r);

end

