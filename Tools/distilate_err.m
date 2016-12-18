function [measures_err, best_err] = distilate_err(Z, param1_name, skinny)

[~, N] = size(Z);

if nargin < 3
    skinny = true;
end

measures_err.param1 = zeros(1, N);
measures_err.psnr = zeros(1, N);
measures_err.fsim = zeros(1,N);
measures_err.ssim = zeros(1,N);
measures_err.msam = zeros(1,N);
measures_err.rel_norm = zeros(1, N);
measures_err.nnz = zeros(1, N);

for j=1:N
    if ~isa(Z{j}, 'MException')
        measures_err.param1(j) = Z{j}.(param1_name);
    else
        if j == 1
            measures_err.param1(j) = eps;
        else
            measures_err.param1(j) = measures_err.param1(j-1) + eps;
        end
    end
    
    for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
        n = i{1};
        measures_err.(n)(j) = get_measure_err(Z{j}, n);
    end
    measures_err.nnz(j) = Z{j}.nnz_e;
end
    
for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
    n = i{1};
    if ~strcmp('msam', n) && ~strcmp('rel_norm', n)
        [best_err.(n).value, best_err.(n).index] = max(measures_err.(n));
    else
        [best_err.(n).value, best_err.(n).index] = min(measures_err.(n));
    end
    best_err.(n).param1 = measures_err.param1(best_err.(n).index);
    if ~skinny && ~isa(Z{best_err.(n).index}, 'MException')
        best_err.(n).L = Z{best_err.(n).index}.L;
        best_err.(n).E = Z{best_err.(n).index}.E;
    end
    best_err.(n).err = Z{best_err.(n).index}.err;
    best_err.(n).iter = Z{best_err.(n).index}.iter;
    best_err.(n).time = Z{best_err.(n).index}.time;
end

end
