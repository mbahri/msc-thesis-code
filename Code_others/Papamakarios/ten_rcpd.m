function [Y, info] = ten_rcpd(X, options)
% [Y, info] = ten_rcpd(X, options)
% Performs Robust CANDECOMP/PARAFAC Decomposition on tensor X.
% INPUTS
%       X           a tensor of any dimensionality
%       options     a struct with options (optional)
% OUTPUTS
%       A           low-rank tensor in Kruskal form
%       E           the sparse component of X
%       info        information about execution
%
% Georgios Papamakarios
% Imperial College London
% May 2014

X = tensor(X);

% options
if nargin < 2
    options = struct;
end
if ~isfield(options, 'return_error') 
    options.return_error = true; 
end
if ~isfield(options, 'return_factors') 
    options.return_factors = false; 
end
if ~isfield(options, 'method') 
    options.method = 'sub'; 
end
options.store_info = nargout >= 2;

switch options.method
    case 'lin'
        [A, E, info] = ten_rcpd_lin(X, options);
    case 'sub'
        [A, E, info] = ten_rcpd_sub(X, options);
    otherwise
        error('Unknown method.');
end

if options.return_factors
    Y.A = A;
else
    Y.A = tensor(A);
end
if options.return_error
    Y.E = E;
end


function [A, E, info] = ten_rcpd_lin(X, options)
% Performs Robust CANDECOMP/PARAFAC Decomposition on tensor X 
% using linearisation.
% INPUTS
%       X           a tensor of any dimensionality
%       options     a struct with options (optional)
% OUTPUTS
%       A           low-rank tensor in Kruskal form
%       E           the sparse component of X
%       info        information about execution

N = ndims(X);
dims = size(X);

% options
if ~isfield(options, 'error_weights') 
    options.error_weights = 1;
end
if ~isfield(options, 'p') 
    options.p = 1; 
end
if ~isfield(options, 'q') 
    options.q = 1; 
end
if ~isfield(options, 'numComp') 
    options.numComp = min(dims); 
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
    U = cell(1, N);
    for i = 1:N
        U{i} = rand(dims(i), options.numComp);
    end
    options.initA = ktensor(U);
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
if ~isfield(options, 'store_info') 
    options.store_info = false; 
end

% initialisation
W = double(options.error_weights);
p = options.p;
q = options.q;
A = options.initA;
E = options.initE;
L = tenzeros(dims);
a = options.alpha;
l = options.lambda;
m = options.initMu;
m_max = options.maxMu;
r = options.rho;
err = inf;
iter = 0;
datanorm = norm(X);
termCrit = options.delta * datanorm;
info.iter = [];
info.err  = [];

% iterating
while err > termCrit && iter < options.maxIter

    % update
    for i = 1:N
        UiT = khatri_rao([A.U(end:-1:i+1) A.U(i-1:-1:1)]);
        n = max(norm(UiT'*UiT, 'fro'), 1.0);
        A.U{i} = sv_shrinkage(A.U{i} - (A.U{i}*UiT' - double(tenmat(X - E + L/m, i))) * UiT / n, a(i)/(m*n), p);
    end
    tenA = tensor(A);
    E = ten_shrinkage(X - tenA + L/m, W*l/m, q);

    % Lagrange multipliers and penalty parameter
    L = L + m * (X - tenA - E);
    m = min(r*m, m_max);

    % calculate error
    err = norm(X - tenA - E);
    iter = iter + 1;

    % print info
    if options.verbose
        fprintf('Iteration %d, error = %g \n', iter, err/datanorm);
    end
    
    % store info
    if options.store_info
        info.iter = [info.iter iter];
        info.err  = [info.err  err/datanorm];
    end

end


function [A, E, info] = ten_rcpd_sub(X, options)
% Performs Robust CANDECOMP/PARAFAC Decomposition on tensor X 
% using substitution.
% INPUTS
%       X           a tensor of any dimensionality
%       options     a struct with options (optional)
% OUTPUTS
%       A           low-rank tensor in Kruskal form
%       E           the sparse component of X
%       info        information about execution

N = ndims(X);
dims = size(X);

% options
if ~isfield(options, 'error_weights') 
    options.error_weights = 1;
end
if ~isfield(options, 'p') 
    options.p = 1; 
end
if ~isfield(options, 'q') 
    options.q = 1; 
end
if ~isfield(options, 'numComp') 
    options.numComp = min(dims); 
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
    U = cell(1, N);
    for i = 1:N
        U{i} = rand(dims(i), options.numComp);
    end
    options.initA = ktensor(U);
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
if ~isfield(options, 'store_info') 
    options.store_info = false; 
end

% initialisation
W = double(options.error_weights);
p = options.p;
q = options.q;
k = options.numComp;
J = cell(1, N);
A = options.initA;
E = options.initE;
L1 = tenzeros(dims);
L2 = cell(1, N);
for i = 1:N
    J{i}  = A.U{i};
    L2{i} = zeros(dims(i), k);
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
info.iter = [];
info.err  = [];

% iterating
while err > termCrit && iter < options.maxIter

    % update
    for i = 1:N
        J{i} = sv_shrinkage(A.U{i} + L2{i}/m, a(i)/m, p);
        UiT = khatri_rao([A.U(end:-1:i+1) A.U(i-1:-1:1)]);
        A.U{i} = (double(tenmat(X - E + L1/m, i)) * UiT + J{i} - L2{i}/m) / (UiT'*UiT + eye(k));
    end
    tenA = tensor(A);
    E = ten_shrinkage(X - tenA + L1/m, W*l/m, q);

    % Lagrange multipliers and penalty parameter
    L1 = L1 + m * (X - tenA - E);
    for i = 1:N
        L2{i} = L2{i} + m * (A.U{i} - J{i});
    end
    m = min(r*m, m_max);

    % calculate error
    err1 = norm(X - tenA - E);
    err2 = zeros(1, N);
    for i = 1:N
        err2(i) = norm(A.U{i} - J{i}, 'fro');
    end
    err = max([err1 err2]);
    iter = iter + 1;

    % print info
    if options.verbose
        fprintf('Iteration %d, error = %g \n', iter, err/datanorm);
    end
    
    % store info
    if options.store_info
        info.iter = [info.iter iter];
        info.err  = [info.err  err/datanorm];
    end

end
