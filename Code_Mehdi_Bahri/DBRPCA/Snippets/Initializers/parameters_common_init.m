function [ vars, params ] = parameters_common_init( X, params, varargin )
%PARAMETERS_COMMON_INIT Default parameters and variable initialization
%   Initialize the variables and parameters shared by all algorithms to
%   common default values. Specificities can be added/overriden in the
%   variants' specific files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

% Dimensions of the input
params.n = size(X, 1);
params.m = size(X, 2);
params.Nobs = size(X, 3);

% Upper bound on the rank is by default the min of the width/height of the
% slices
% params.r = min(params.n, params.m);
params = set_if_unset(params, {'r'}, {min(params.n, params.m)});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters

params.MONOTONE = false;
params.PARALLEL = false;
params.TIME = false;
params.MEAN = false;
params.MMX = true;

% params.lambda = 1 / sqrt(params.Nobs * max(params.n, params.m));
% 
% % params.lambda = 1 / sqrt(max([params.n, params.m]));
% params.tol = 1e-7;
% params.MAXITER = 300;
% params.rho = 1.2;
% 
% params.alpha_c = 1;
% params.alpha_r = 1;
% params.alpha_t = 1e-1;

params = set_if_unset(params, {...
        'lambda', ...
        'tol', ...
        'MAXITER', ...
        'rho', ...
        'alpha_c', ...
        'alpha_r', ...
        'alpha_t'...
    }, ...
    {...
        1 / sqrt(params.Nobs * max(params.n, params.m)), ...
        1e-7, ...
        150, ...
        1.2, ...
        1, ...
        1, ...
        1e-1...
    });

params.mu = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the optional parameters

params = parse_input_parameters(varargin, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable definitions
% Input
vars.X = X;
% Sparse error terms
vars.E = zeros(params.n, params.m, params.Nobs);
% Lagrange multipliers
vars.Y = zeros(params.n, params.m, params.Nobs);

% Robust mean
if params.MEAN
    vars.M = zeros(params.n, params.m, params.Nobs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update functions definition
params = set_update_functions_l2(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization with SVD
vars = init_via_svd(vars, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Penalty parameter
% Initialization of mu that is the same for all algorithms
if params.mu < 0
    params.mu = init_mu_with_norms(vars.X, 1.25, params);
end
% Default upper bound on the penalty
params.mu_bar = inf;
% params.mu_bar = 1e9*params.mu;
end

function [st] = set_if_unset(st, names, values)

for i=1:length(names)
    name = names{i};
    value = values{i};

    if ~isfield(st, name)
        st.(name) = value;
    end
    
end

end