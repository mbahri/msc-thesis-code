function [ params ] = set_update_functions_l2( params )
%SET_UPDATE_FUNCTIONS_L2 Sets the update functions depending on input params
%
% Mehdi Bahri - Imperial College London
% July, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stopping criterion
params.stopping_criterion = @(vars, params) (...
    stop_max_error(vars, params) ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behaviour that depends on MMX
if params.MMX
    params.comp_Uc = @(vars, params) (...
        comp_Uc_mmx_L2(vars, params)...
    );
    params.comp_Ur = @(vars, params) (...
        comp_Ur_mmx_L2(vars, params)...
    );
    params.update_E = @(vars, params) ( ...
        update_E_mmx(vars, params) ...
    );
    params.update_TY = @(vars, params) (...
        update_TY_mmx_regT_L2(vars, params) ...
    );
else
    params.comp_Uc = @(vars, params) (...
        comp_Uc_nommx_L2(vars, params)...
    );
    params.comp_Ur = @(vars, params) (...
        comp_Ur_nommx_L2(vars, params)...
    );
    params.update_E = @(vars, params) ( ...
        update_E_nommx(vars, params) ...
    );
    params.update_TY = @(vars, params) (...
        update_TY_nommx_regT_L2(vars, params) ...
    );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behaviour that depends on the mean
if params.MEAN
    params.visualize = @(vars, params) (...
        visualize_with_mean(vars, params)...
    );
else
    params.visualize = @(vars, params) (...
        visualize_no_mean(vars, params)...
    );
end

end
