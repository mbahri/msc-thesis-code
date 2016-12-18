function [results] = benchmark_mod(algorithm, data, params, savedata, save_idx, file)
% [results] = benchmark(algorithm, data, params, savedata, file)
% Runs a particular method for a set of parameters and stores the results.
% INPUTS
%       algorithm   a function handle to the algorithm to run
%       data        a struct containing the data and the low-rank component
%       params      a struct with lists of parameters to run the algorithm for
%       savedata    if true, also stores the recovered low-rank component (optional)
%       file        if specified, saves the result to this file (optional)
% OUTPUTS
%       results     a struct containing the results for all parameters
%
% Georgios Papamakarios
% Imperial College London
% June 2014
%
% Modified by Mehdi Bahri

if nargin < 6 
    file = [];
end

if nargin < 5
    save_idx = 1:size(data.X, 3);
end

if nargin < 4 
    savedata = true;
end

if ~isfield(params, 'numComp')
    params.numComp = {-1};
end 
if ~iscell(params.numComp)
    params.numComp = num2cell(params.numComp);
end

len_l = length(params.lambda);
len_p = length(params.p);
len_r = length(params.numComp);
   
results.err  = zeros(len_l, len_p, len_r);
results.time = zeros(len_l, len_p, len_r);
results.info = cell(len_l, len_p, len_r);
if savedata
    results.A = cell(len_l, len_p, len_r);
end

options.return_error = false;
if isfield(params, 'method')
    options.method = params.method;
end 
if isfield(params, 'alpha')
    options.alpha = params.alpha;
end 
if isfield(params, 'error_weights')
    options.error_weights = params.error_weights;
end 

for i = 1:len_l
    for j = 1:len_p
        for k = 1:len_r
        
            options.lambda = params.lambda(i);
            options.p = params.p(j);
            options.q = params.p(j);
            options.numComp = params.numComp{k};

            try
            
                % run the algorithm
                tic
                [Y, results.info{i,j,k}] = algorithm(data.X, options);
                results.time(i,j,k) = toc;

                if savedata
                    results.A{i,j,k} = Y.A(:,:,save_idx);
                end
                if isfield(options, 'error_weights')
                    idx = ~logical(options.error_weights);
                    dA = double(data.A);
                    YA = double(Y.A);
                    results.err(i,j,k) = ew_norm(dA(idx) - YA(idx), 2) / ew_norm(dA(idx), 2);
                else  
                    results.err(i,j,k) = ew_norm(data.A - Y.A, 2) / ew_norm(data.A, 2);
                end
            
            catch exception
                
                fprintf('ERROR (lambda = %g, p = %g, numComp = [ ', options.lambda, options.p);
                fprintf('%g ', options.numComp);
                fprintf(']): %s \n', exception.message);
                
                results.info{i,j,k} = nan;
                results.time(i,j,k) = nan;
                if savedata
                    results.A{i,j,k} = nan;
                end
                results.err(i,j,k) = nan;
            
            end
        end
    end
end

% save to file (optionally)
if ~isempty(file)
    save(file, 'method', 'params', 'results', '-v7.3');
end
