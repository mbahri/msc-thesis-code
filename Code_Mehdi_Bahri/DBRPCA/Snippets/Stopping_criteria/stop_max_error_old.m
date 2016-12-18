function [ emax ] = stop_max_error( vars, params )
%STOP_MAX_ERROR Computes the maximum reconstruction error of the slices
%
% Mehdi Bahri - Imperial College London
% April, 2016

X = vars.X;
Xt = vars.Xt;   % Xt is X - E or X - E - M as needed
T = vars.T;
Uc = vars.Uc;
Ur = vars.Ur;

emax = -inf;
for k=1:params.Nobs
    emax = get_max_error(X(:,:,k), Xt(:,:,k), T(:,:,k), ...
        Uc, Ur, emax);
end

end

function [emax] = get_max_error(X, Xt, T, Uc, Ur, currmax)

Error = norm(Xt - Uc*T*Ur', 'fro') /...
    norm(X, 'fro');

if Error > currmax
   emax = Error;
else
   emax = currmax;
end

end