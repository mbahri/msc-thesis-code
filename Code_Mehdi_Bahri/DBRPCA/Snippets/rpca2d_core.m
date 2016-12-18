function [ Data, Info ] = rpca2d_core( vars, params)
%RPCA2D_CORE Core ADMM algorithm for NO 2D RPCA
%
% Mehdi Bahri - Imperial College London
% April, 2016


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convergence
converged = false;
niter = 0;
EE = inf;
err = [];   % Store the errors for plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm
params

while ~converged && niter < params.MAXITER
    if params.TIME > 0
        tic
    end
    niter = niter + 1;
    
    % ------ Update En first ------
    vars = params.update_E(vars, params);
    
    % ------ Robust mean if enabled ------
    if params.MEAN
        vars.M = sum(vars.X - vars.E, 3) / params.Nobs;
        vars.M = repmat(vars.M, 1, 1, params.Nobs);
        
        % "X tilde" is X - E - M
        vars.Xt = vars.X - vars.E - vars.M;
    else
        % "X tilde" is X - E
        vars.Xt = vars.X - vars.E;
    end
    
    % ------ Start with Uc and Ur ------
    vars.S = params.mu * vars.Xt + vars.Y;
    vars = params.comp_Uc(vars, params);
    vars = params.comp_Ur(vars, params);

    % ------ Solve for each T ------
    vars = params.update_TY(vars, params);
    
    % ------ Update mu ------
    params.mu = min(params.mu_bar, params.rho * params.mu);
    
    % ------ Test for convergence and monotonicity ------
    EEp = EE;
    EE = params.stopping_criterion(vars, params);
    if params.MONOTONE && EEp < EE
        fprintf('Error increasing - stopping now\n');
        converged = true;
    end
    if EE < params.tol
        converged = true;
    end
    
    % Update the errors
    err = [err EE];
    
    % ------ Some profiling information ------
    if params.TIME > 0
        itertime = toc;
    else
        itertime = nan;
    end
    fprintf('[%03d] mu = %f max error = %g | rk(Uc, 1e-3) = %d rk(Ur, 1e-3) = %d | %fs\n', ...
            niter, params.mu, EE, estim_rank(vars.Uc, 1e-3), estim_rank(vars.Ur, 1e-3), itertime);
    params.visualize(vars, params);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Returned variables
Data = vars;

Info.niter = niter;
Info.err = err;

end