function [ emax ] = stop_max_error( vars, params )
%STOP_MAX_ERROR Computes the maximum reconstruction error of the slices
%
% Mehdi Bahri - Imperial College London
% April, 2016

X = vars.X;
Xt = vars.Xt;   % Xt is X - E or X - E - M as needed

% emax = -inf;
e_slice = -inf;
e_base = -inf;
e_core = -inf;

if isfield(vars, 'A')
    Uc = vars.A;
    Ur = vars.B;
else
    Uc = vars.Uc;
    Ur = vars.Ur;
end

if isfield(vars, 'K')
    T = vars.K;
else
    T = vars.T;
end

for k=1:params.Nobs
    e_slice = get_max_error(X(:,:,k), Xt(:,:,k), T(:,:,k), ...
        Uc, Ur, e_slice);
end

if isfield(vars, 'A')
    e_base = max(matrix_relative_error(vars.A, vars.Uc), matrix_relative_error(vars.B, vars.Ur));
end

if isfield(vars, 'K')
    for k=1:params.Nobs
        e_core = max(e_core, matrix_relative_error(vars.K(:,:,k), vars.T(:,:,k)));
    end
end

emax = max([e_slice, e_base, e_core]);

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