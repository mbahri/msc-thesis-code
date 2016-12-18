function [measures_mean, best_mean] = distilate_mean(Z, param1_name, skinny)

[~, N] = size(Z);

if nargin < 3
    skinny = true;
end

measures_mean.param1 = zeros(1, N);
measures_mean.psnr = zeros(1, N);
measures_mean.fsim = zeros(1,N);
measures_mean.ssim = zeros(1,N);
measures_mean.msam = zeros(1,N);
measures_mean.rel_norm = zeros(1, N);
measures_mean.nnz = zeros(1, N);

for j=1:N
    if ~isa(Z{j}, 'MException')
        measures_mean.param1(j) = Z{j}.(param1_name);
    else
        if j == 1
            measures_mean.param1(j) = eps;
        else
            measures_mean.param1(j) = measures_mean.param1(j-1) + eps;
        end
    end
    
    for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
        n = i{1};
        measures_mean.(n)(j) = get_measure_err(Z{j}, n);
    end
    measures_mean.nnz(j) = Z{j}.nnz_e;
end
    
for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
    n = i{1};
    if ~strcmp('msam', n) && ~strcmp('rel_norm', n)
        [best_mean.(n).value, best_mean.(n).index] = max(measures_mean.(n));
    else
        [best_mean.(n).value, best_mean.(n).index] = min(measures_mean.(n));
    end
    best_mean.(n).param1 = measures_mean.param1(best_mean.(n).index);
    if ~skinny && ~isa(Z{best_mean.(n).index}, 'MException')
        best_mean.(n).L = Z{best_mean.(n).index}.L;
        best_mean.(n).E = Z{best_mean.(n).index}.E;
    end
    best_mean.(n).err = Z{best_mean.(n).index}.err;
    best_mean.(n).iter = Z{best_mean.(n).index}.iter;
    best_mean.(n).time = Z{best_mean.(n).index}.time;
end

end
