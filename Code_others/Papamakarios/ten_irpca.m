function [Y, info] = ten_irpca(X, options)
% [Y, info] = ten_irpca(X, options)
% Performs Inductive Robust Principal Component Analysis on tensor X.
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
        [P, E, info] = ten_irpca_lin(X, options);
    case 'sub'
        [P, E, info] = ten_irpca_sub(X, options);
    otherwise
        error('Unknown method.');
end

if options.return_factors
    Y.P = P;
else
    Y.A = ttm(X, P);
end
if options.return_error
    Y.E = E;
end


function [P, E, info] = ten_irpca_lin(X, options)
% Performs Inductive Robust Principal Component Analysis on tensor X
% using linearisation.
% INPUTS
%       X           a tensor of any dimensionality
%       options     a struct with options
% OUTPUTS
%       P           cell array with the projection matrices for each mode
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
    options.numComp = floor(0.5 * dims); 
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
if ~isfield(options, 'initP')  
    options.initP = cell(1, N);
    for i = 1:N
        Ui = nvecs(X, i, options.numComp(i));
        options.initP{i} = Ui * Ui';
    end
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
P = options.initP;
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

    % updates
    for i = 1:N
        XPi = double(tenmat(ttm(X, P, -i), i));
        n = max(norm(XPi*XPi', 'fro'), 1.0);
        P{i} = sv_shrinkage(P{i} - (P{i}*XPi - double(tenmat(X - E + L/m, i))) * XPi' / n, a(i)/(m*n), p);
    end
    E = ten_shrinkage(X - ttm(X,P) + L/m, W*l/m, q);

    % Lagrange multipliers and penalty parameter
    L = L + m * (X - ttm(X,P) - E);
    m = min(r*m, m_max);

    % calculate error
    err = norm(X - ttm(X,P) - E);
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


function [P, E, info] = ten_irpca_sub(X, options)
% Performs Inductive Robust Principal Component Analysis on tensor X 
% using substitution.
% INPUTS
%       X           a tensor of any dimensionality
%       options     a struct with options
% OUTPUTS
%       P           cell array with the projection matrices for each mode
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
    options.numComp = floor(0.5 * dims); 
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
if ~isfield(options, 'initP')  
    options.initP = cell(1, N);
    for i = 1:N
        Ui = nvecs(X, i, options.numComp(i));
        options.initP{i} = Ui * Ui';
    end
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
J = cell(1, N);
P = options.initP;
E = options.initE;
L1 = tenzeros(dims);
L2 = cell(1, N);
for i = 1:N
    J{i}  = P{i};
    L2{i} = zeros(dims(i));
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

    % updates
    for i = 1:N
        J{i} = sv_shrinkage(P{i} + L2{i}/m, a(i)/m, p);
        XPi = double(tenmat(ttm(X, P, -i), i));
        P{i} = (double(tenmat(X - E + L1/m, i)) * XPi' + J{i} - L2{i}/m) / (XPi*XPi' + eye(size(X,i)));
    end
    E = ten_shrinkage(X - ttm(X,P) + L1/m, W*l/m, q);

    % Lagrange multipliers and penalty parameter
    L1 = L1 + m * (X - ttm(X,P) - E);
    for i = 1:N
        L2{i} = L2{i} + m * (P{i} - J{i});
    end
    m = min(r*m, m_max);

    % calculate error
    err1 = norm(X - ttm(X,P) - E);
    err2 = zeros(1, N);
    for i = 1:N
        err2(i) = norm(P{i} - J{i}, 'fro');
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
