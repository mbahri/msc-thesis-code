function [ MSE ] = MSE( Ori, Dist )
%MSE Mean squared error between two samples
%   Mehdi Bahri - Imperial College London, July 2016

[M N] = size(Ori);
error = Ori - Dist;
MSE = sum(sum(error .* error)) / (M * N);


end

