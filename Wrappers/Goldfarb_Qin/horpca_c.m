function [ Data, Info ] = horpca_c( X, r, mu, isMV )
%HORPCA_C Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    isMV = false;
end

%%%%%%%% DATA
data.X = tensor(X);
data.T = tensor(X);
data.b = X(:);

if isMV
   data.linInd = find(X ~= 0);
   data.b = data.b(data.linInd);
else
    data.linInd = (1:numel(X))';
end

%%%%%%%% PARAMETERS
params.X0 = tenzeros( size(data.T) );
N = ndims(data.T);
params.V0 = cell( 1, N );
for i = 1:N
    params.V0{i} = tenzeros( size(data.T) );
end
params.E0 = tenzeros( size(data.T) );
params.mu0 = 1/(N+1);

mu1fac = 10; % Was 5 in the script but Appendix E says 10
params.mu1fac = mu1fac;

% For non-convex, we choose mu = 1 like in Appendix E
params.mu1 = mu;  %mu1fac*std(T(:));
params.mu2 = params.mu1;
params.mu_min = 1e-4;   % needed?
params.mu_max = 1e2;    % needed?

params.max_iter = 1000;
params.opt_tol = 1e-3;  %1e-5 for syn data analysis plots, 1e-3 for others
params.eta = 1/(N+1);   % needed?
params.IsTC = false;    % needed?

% Upper bound on the Tucker rank
params.k = r;
params.verbose = true;
params.use_cont = true; % Needed outside of FISTA?

%%%%%%%%%% for PROPACK %%%%%%%%%%%%
% declare global var 'sv'
global sv;
global tmode;
global use_propack;
global curr_mu;
sv =  ceil(min(size(data.T)) * 0.1) * ones( 1, N );
use_propack = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% ACTUALLY RUN THE ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%

warning('off','PROPACK:NotUsingMex');

if isMV
    results = tensor_rpca_tc_adal_ncx( data, params );
else
    results = tensor_rpca_adal_ncx( data, params );
end

Data.L = double(results.X);
Data.E = double(results.E);
Data.U = results.U;

Info.err = results.err;
Info.iter = 1:results.iter;

end

