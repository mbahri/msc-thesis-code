function [ Data, Info ] = rpca2d_l2( X, varargin)
%RPCA2D_L2 Solves NO2DRPCA with a L2 penalty on T and Frobenius norm
% min (1/2)*(alpha_c*||Uc||_F^2 + alpha_c*||Ur||_F^2 + 
%                   alpha_r*sum(||Tn||_F^2) + sum(labmda_n*||En||_1)
% s.t for all n Xn = UcTnUr^T + En
%
% Mehdi Bahri - Imperial College London
% April, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults for this algorithm
params = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shared behaviour and default parameter initalization
[vars, params] = parameters_common_init( X, params, varargin{:} );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm and specific output values
[Out, Info] = rpca2d_core(vars, params);

Data.Uc = Out.Uc;
Data.Ur = Out.Ur;
Data.T = Out.T;
Data.E = Out.E;

if params.MEAN
    Data.M = Out.M(:,:,1);
    Data.L = wm_make_L(Out.Uc, Out.Ur, Out.T) + Out.M;
else
    Data.L = wm_make_L(Out.Uc, Out.Ur, Out.T);
end

end