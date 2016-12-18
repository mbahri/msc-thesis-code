function [results] = benchmark_all_ten(data, lambda, p, nrank, rank, savedata, file, alpha, error_weights)
% [results] = benchmark_all_ten(data, lambda, p, nrank, rank, savedata, file, alpha, error_weights)
% Benchmarks all tensor algorithms on given data and for a given set of
% parameters and stores the result.
% INPUTS
%       data            a struct containing the data tensor and its low-rank component
%       lambda          an array of lambda's to run for
%       p               an array of p (= q) to run for
%       nrank           an array of n-ranks to run for
%       rank            an array of ranks to run for 
%       savedata        if true, also stores the recovered low-rank component (optional)
%       file            if specified, saves the result to this file (optional)
%       alpha           the alpha values to use for the algorithms (optional)
%       error_weights   the error weights to use for the algorithms (optional)
% OUTPUTS
%       results         a struct containing the results for all parameters
%
% Georgios Papamakarios
% Imperial College London
% June 2014

if nargin < 7
    file = [];
end
if nargin < 6
    savedata = true;
end

params.lambda = lambda;
params.p = p;
if nargin > 8
    params.error_weights = error_weights;
end
if nargin > 7
    params.alpha = alpha;
end

% RPCA
fprintf('RPCA \n');
if isfield(params, 'numComp'), params = rmfield(params, 'numComp'); end
results.rpca = benchmark(@ten_rpca, data, params, savedata);

% BRPCA
fprintf('BRPCA \n');
params.numComp = nrank;
results.brpca = benchmark(@ten_brpca, data, params, savedata);

% IRPCA sub
fprintf('IRPCA (sub) \n');
params.method = 'sub';
params.numComp = nrank;
results.irpca.sub = benchmark(@ten_irpca, data, params, savedata);

% IRPCA lin
fprintf('IRPCA (lin) \n');
params.method = 'lin';
params.numComp = nrank;
results.irpca.lin = benchmark(@ten_irpca, data, params, savedata);

% ORPCA
fprintf('ORPCA \n');
params.numComp = nrank;
results.orpca = benchmark(@ten_orpca, data, params, savedata);

% RCPD sub
fprintf('RCPD (sub) \n');
params.method = 'sub';
params.numComp = rank;
results.rcpd.sub = benchmark(@ten_rcpd, data, params, savedata);

% RCPD lin
fprintf('RCPD (lin) \n');
params.method = 'lin';
params.numComp = rank;
results.rcpd.lin = benchmark(@ten_rcpd, data, params, savedata);

% save to file (optionally)
if ~isempty(file)
    save(file, 'results', '-v7.3');
end
