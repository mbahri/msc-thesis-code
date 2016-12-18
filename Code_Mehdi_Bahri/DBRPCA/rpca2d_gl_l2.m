function [ Data, Info ] = rpca2d_gl_l2( X, varargin)
%RPCA2D_GL_L2 Solves NO2DRPCA with a L2 penalty on T and Group Lasso
% min (1/2)*(alpha_c*||Uc||_F^2 + alpha_c*||Ur||_F^2 + 
%                   alpha_r*sum(||Tn||_F^2) + sum(labmda_n*||En||_1)
% s.t for all n Xn = UcTnUr^T + En
%
% Mehdi Bahri - Imperial College London
% July, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults for this algorithm
params = struct();
params.alpha_c = 1e-3;
params.alpha_r = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shared behaviour and default parameter initalization
[vars, params] = parameters_common_init( X, params, varargin{:} );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update functions overrides
params.comp_Uc = @(vars, params) (...
    comp_Uc_mmx_group_lasso(vars, params)...
);
params.comp_Ur = @(vars, params) (...
    comp_Ur_mmx_group_lasso(vars, params)...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional variables and overrides
vars.A = vars.Uc;
vars.B = vars.Ur;

vars.Y_c = zeros(params.n, params.r);
vars.Y_r = zeros(params.m, params.r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional penalty parameters and overrides

% params.mu_c = 1.25 / norm(vars.Uc, 'fro');
% params.mu_r = 1.25 / norm(vars.Ur, 'fro');
vars.mu_c = 1e-3;
vars.mu_r = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm and specific output values
[Out, Info] = rpca2d_core(vars, params);

Data.Uc = Out.A;
Data.Ur = Out.B;
Data.T = Out.T;
Data.E = Out.E;

if params.MEAN
    Data.M = Out.M(:,:,1);
    Data.L = wm_make_L(Out.Uc, Out.Ur, Out.T) + Out.M;
else
    Data.L = wm_make_L(Out.Uc, Out.Ur, Out.T);
end 
%     % Test updating Ur with the previous iteration of Uc (i.e before Uc is
%     % updated)
%     Uc = vars.Uc;
%     B = vars.B;
%     [vars.A, vars.Uc, vars.Y_c, params.mu_c] = comp_Uc_mmx_group_lasso(S, vars.T, vars.Ur, vars.A, vars.Y_c, params);
%     [vars.B, vars.Ur, vars.Y_r, params.mu_r] = comp_Ur_mmx_group_lasso(S, vars.T, Uc, B, vars.Y_r, params);
    
%     % Test enforcing the column sparse Uc, Ur
%     [vars.Uc, vars.A, vars.Y_c, params.mu_c] = comp_Uc_mmx_group_lasso(S, vars.T, vars.B, vars.Uc, vars.Y_c, params);
%     [vars.Ur, vars.B, vars.Y_r, params.mu_r] = comp_Ur_mmx_group_lasso(S, vars.T, vars.A, vars.Ur, vars.Y_r, params);

end