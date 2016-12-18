function [ L, E ] = wrapper_mehdi( fun, X, r, lambda, alpha_t )
%WRAPPER_MEHDI Calls Mehdi Bahri's code and returns L and E.
%
% Mehdi Bahri - Imperial College London
% July, 2016

if nargin < 4
    [Uc, Ur, T, E] = fun(X, r);
elseif nargin < 5
    [Uc, Ur, T, E] = fun(X, r, 'lambda', lambda);
else
    [Uc, Ur, T, E] = fun(X, r, 'lambda', lambda, 'alpha_t', alpha_t);
end
    
L = wm_make_L(Uc, Ur, T);

% L = mmx('mult', Uc, T);
% L = mmx('mult', L, Ur, 'nt');

end

