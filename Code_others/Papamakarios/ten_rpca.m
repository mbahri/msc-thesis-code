function [Y, info] = ten_rpca(X, options)
% [Y, info] = ten_rpca(X, options)
% Performs Robust Principal Component Analysis on tensor X.
% INPUTS
%       X           a tensor of any dimensionality
%       options     a struct with options (optional)
% OUTPUTS
%       Y           a struct containing the result
%       info        information about execution
%
% Georgios Papamakarios
% Imperial College London
% May 2014

X = tensor(X);
N = ndims(X);
dims = size(X);

% options
if nargin < 2
    options = struct;
end
if ~isfield(options, 'error_weights') 
    options.error_weights = 1;
end
if ~isfield(options, 'return_error') 
    options.return_error = true; 
end
if ~isfield(options, 'p') 
    options.p = 1; 
end
if ~isfield(options, 'q') 
    options.q = 1; 
end
if ~isfield(options, 'alpha') 
    options.alpha = ones(1,N) / N; 
end
if ~isfield(options, 'lambda') 
    options.lambda = 1 / sqrt(max(dims));
end
if ~isfield(options, 'delta')  
    options.delta = 1.0e-7; 
end
if ~isfield(options, 'initMu')  
    options.initMu = 1.0e-3; 
end
if ~isfield(options, 'maxMu')  
    options.maxMu = 1.0e+9; 
end
if ~isfield(options, 'rho')  
    options.rho = 1.2; 
end
if ~isfield(options, 'initA')  
    options.initA = tenzeros(dims); 
end
if ~isfield(options, 'initE')  
    options.initE = tenzeros(dims); 
end
if ~isfield(options, 'verbose')  
    options.verbose = false; 
end
if ~isfield(options, 'maxIter')  
    options.maxIter = inf; 
end

% initialisation
W = double(options.error_weights);
p = options.p;
q = options.q;
J = cell(1, N);
A = options.initA;
E = options.initE;
L1 = tenzeros(dims);
L2 = cell(1, N);
for i = 1:N
    J{i}  = tenmat(A, i);
    L2{i} = tenmat(tenzeros(dims), i);
end    
a = options.alpha;
l = options.lambda;
m = options.initMu;
m_max = options.maxMu;
r = options.rho;
err = inf;
iter = 0;
datanorm = norm(X);
termCrit = options.delta * datanorm;
store_info = nargout >= 2;
info.iter = [];
info.err  = [];

% iterating
while err > termCrit && iter < options.maxIter

    % updates
    sumJ = tenzeros(dims);
    for i = 1:N
        J{i}(:,:) = sv_shrinkage(tenmat(A,i) + L2{i}*(1/m), a(i)/m, p);
        sumJ = sumJ + tensor(J{i} - L2{i}*(1/m));
    end
    A = (X - E + L1/m + sumJ) / (N+1);
    E = ten_shrinkage(X - A + L1/m, W*l/m, q);

    % Lagrange multipliers and penalty parameter
    L1 = L1 + m * (X - A - E);
    for i = 1:N
        L2{i} = L2{i} + m * (tenmat(A,i) - J{i});
    end
    m = min(r*m, m_max);

    % calculate error
    err1 = norm(X - A - E);
    err2 = zeros(1, N);
    for i = 1:N
        err2(i) = norm(tenmat(A,i) - J{i});
    end
    err = max([err1 err2]);
    iter = iter + 1;

    % print info
    if options.verbose
        fprintf('Iteration %d, error = %g \n', iter, err/datanorm);
    end
    
    % store info
    if store_info
        info.iter = [info.iter iter];
        info.err  = [info.err  err/datanorm];
    end

end

Y.A = A;
if options.return_error
    Y.E = E;
end
