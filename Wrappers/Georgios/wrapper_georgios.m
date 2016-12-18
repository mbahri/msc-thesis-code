function [ L, E ] = wrapper_georgios( fun, X, opts )
%WRAPPER_GEORGIOS Calls Georgios Papamakario's functions and returns L and
%E.
%
% Mehdi Bahri - Imperial College London
% July, 2016

[Y, ~] = fun(X, opts);
L = double(Y.A);
E = double(Y.E);

end

