function [data] = create_synthetic_data(opt)
% [data] = create_synthetic_data(opt)
% Creates a random triplet of matrices or tensors, composed of a low-rank
% component, a sparse component and their superposition.
% INPUTS
%       opt    a struct of options
% OUTPUTS
%       data   a struct containing the data
%
% Georgios Papamakarios
% Imperial College London
% June 2014

if ~isfield(opt, 'return_error')
    opt.return_error = true;
end

switch opt.type
    
    case 'matrix'
        
        % low-rank component
        rnk = floor(opt.rank_perc * min(opt.dims));
        U = randn(opt.dims(1), rnk);
        V = randn(rnk, opt.dims(2));
        A = U * V;
        A = A / std(A(:), 1) * opt.data_magn;
        
        % sparse component
        E = opt.err_magn * (2 * rand(opt.dims) - 1);
        idx = rand(opt.dims) > opt.err_perc;
        E(idx) = 0;
        
    case 'cp_tensor'
        
        % low-rank component
        N = length(opt.dims);
        rnk = floor(opt.rank_perc * min(opt.dims));
        U = cell(1, N);
        for i = 1:N
            U{i} = randn(opt.dims(i), rnk);
        end
        A = full(ktensor(U));
        A = A / std(A(:), 1) * opt.data_magn;
        
        % sparse component
        E = opt.err_magn * (2 * rand(opt.dims) - 1);
        idx = rand(opt.dims) > opt.err_perc;
        E(idx) = 0;
        E = tensor(E);
        
    case 'tucker_tensor'
        
        % low-rank component
        N = length(opt.dims);
        rnk = floor(opt.rank_perc .* opt.dims);
        U = cell(1, N);
        for i = 1:N
            U{i} = randn(opt.dims(i), rnk(i));
        end
        V = tensor(randn(rnk));
        A = ttm(V, U);
        A = A / std(A(:), 1) * opt.data_magn;
        
        % sparse component
        E = opt.err_magn * (2 * rand(opt.dims) - 1);
        idx = rand(opt.dims) > opt.err_perc;
        E(idx) = 0;
        E = tensor(E);
        
    otherwise
        
        error('Unknown data type.');
end

X = A + E;

data.X = X;
data.A = A;
if opt.return_error
    data.E = E;
end
