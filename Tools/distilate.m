function [measures, best] = distilate(Z, param1_name, skinny, ref_fsimc)

[~, N] = size(Z);

if nargin < 3
    skinny = true;
end

measures.param1 = zeros(1, N);
measures.psnr = zeros(1, N);
measures.fsim = zeros(1,N);
measures.ssim = zeros(1,N);
measures.msam = zeros(1,N);
measures.rel_norm = zeros(1, N);

if nargin < 4
    comp_fsimc = false;
else
    O = rescale(ref_fsimc, 0, 255);
    comp_fsimc = true;
end

for j=1:N
    if ~isa(Z{j}, 'MException')
        measures.param1(j) = Z{j}.(param1_name);
    else
        if j == 1
            measures.param1(j) = eps;
        else
            measures.param1(j) = measures.param1(j-1) + eps;
        end
    end
    
    for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
        n = i{1};
        measures.(n)(j) = get_measure(Z{j}, n);
    end
    if comp_fsimc
        if isa(Z{j}, 'MException')
            measures.fsimc(j) = 0;
        else
            [~, measures.fsimc(j)] = FeatureSIM(O, rescale(Z{j}.L, 0, 255));
        end
    end
end

if comp_fsimc
    [best.fsimc.value, best.fsimc.index] = max(measures.fsimc);
    best.fsimc.param1 = measures.param1(best.fsimc.index);
    
    if ~skinny && ~isa(Z{best.fsimc.index}, 'MException')
        best.fsimc.L = Z{best.fsimc.index}.L;
        best.fsimc.E = Z{best.fsimc.index}.E;
    end
%     best.fsimc.err = Z{best.fsimc.index}.err;
    best.fsimc.iter = Z{best.fsimc.index}.iter;
%     best.fsimc.time = Z{best.fsimc.index}.time;
end

for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
    n = i{1};
    if ~strcmp('msam', n) && ~strcmp('rel_norm', n)
        [best.(n).value, best.(n).index] = max(measures.(n));
    else
        [best.(n).value, best.(n).index] = min(measures.(n));
    end
    best.(n).param1 = measures.param1(best.(n).index);
    if ~skinny && ~isa(Z{best.(n).index}, 'MException')
        best.(n).L = Z{best.(n).index}.L;
        best.(n).E = Z{best.(n).index}.E;
    end
%     best.(n).err = Z{best.(n).index}.err;
    best.(n).iter = Z{best.(n).index}.iter;
%     best.(n).time = Z{best.(n).index}.time;
end

end
