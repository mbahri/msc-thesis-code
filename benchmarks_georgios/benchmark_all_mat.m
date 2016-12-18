function [results] = benchmark_all_mat(data, lambda, p, rank, savedata, file)
% [results] = benchmark_all_mat(data, lambda, p, rank, savedata, file)
% Benchmarks all matrix algorithms on given data and for a given set of
% parameters and stores the result.
% INPUTS
%       data        a struct containing the data matrix and its low-rank component
%       lambda      an array of lambda's to run for
%       p           an array of p (= q) to run for
%       rank        an array of ranks to run for 
%       savedata    if true, also stores the recovered low-rank component (optional)
%       file        if specified, saves the result to this file (optional)
% OUTPUTS
%       results     a struct containing the results for all parameters
%
% Georgios Papamakarios
% Imperial College London
% June 2014

if nargin < 6
    file = [];
end
if nargin < 5
    savedata = true;
end

params.lambda = lambda;
params.p = p;

% % CPCA
% fprintf('CPCA \n');
% params.numComp = rank;
% results.cpca = benchmark(@cpca, data, params, savedata);

% RPCA alm
fprintf('RPCA (alm) \n');
params.method = 'alm';
if isfield(params, 'numComp'), params = rmfield(params, 'numComp'); end
results.rpca.alm = benchmark(@mat_rpca, data, params, savedata);

% % RPCA apg
% fprintf('RPCA (apg) \n');
% params.method = 'apg';
% if isfield(params, 'numComp'), params = rmfield(params, 'numComp'); end
% results.rpca.apg = benchmark(@rpca, data, params, savedata);
% 
% % BRPCA
% fprintf('BRPCA \n');
% params.numComp = rank;
% results.brpca = benchmark(@brpca, data, params, savedata);
% 
% % IRPCA sub
% fprintf('IRPCA (sub) \n');
% params.method = 'sub';
% if isfield(params, 'numComp'), params = rmfield(params, 'numComp'); end
% results.irpca.sub = benchmark(@irpca, data, params, savedata);
% 
% % IRPCA lin
% fprintf('IRPCA (lin) \n');
% params.method = 'lin';
% if isfield(params, 'numComp'), params = rmfield(params, 'numComp'); end
% results.irpca.lin = benchmark(@irpca, data, params, savedata);
% 
% % ORPCA
% fprintf('ORPCA \n');
% params.numComp = rank;
% results.orpca = benchmark(@orpca, data, params, savedata);
% 
% % ROSL
% fprintf('ROSL \n');
% params.numComp = rank;
% results.rosl = benchmark(@rosl, data, params, savedata);

% save to file (optionally)
if ~isempty(file)
    save(file, 'results', '-v7.3');
end
