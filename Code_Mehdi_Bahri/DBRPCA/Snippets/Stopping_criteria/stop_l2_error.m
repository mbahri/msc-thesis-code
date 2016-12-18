function [ err ] = stop_l2_error( vars )
%STOP_L2_ERROR Computes the tensor relative reconstruction error
%
% Mehdi Bahri - Imperial College London
% June, 2016

err = norm(tensor(vars.Xt) - ttm(tensor(vars.T), ...
    {vars.Uc, vars.Ur}, [1 2])) / norm(tensor(vars.X));

end

