function [vars] = update_E_mmx(vars, params) %Y, T, Uc, Ur
%UPDATE_E_MMX Solves for E by soft-thresholding, MMX
%
% Mehdi Bahri - Imperial College London
% May, 2016

A = mmx('mult', vars.Uc, vars.T);
A = mmx('mult', A, vars.Ur, 'nt');

% If we estimate the mean we must substract it from X
if params.MEAN
    X = vars.X - vars.M;
else
    X = vars.X;
end

temp_T = X - A + (1/params.mu)*vars.Y;

vars.E = soft_shrinkage(temp_T, params.lambda/params.mu);

end