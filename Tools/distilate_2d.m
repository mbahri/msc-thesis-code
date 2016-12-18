function [measures, best] = distilate_2d(Z, param1_name, param2_name, ref_fsimc)

if nargin < 4
    comp_fsimc = false;
else
    O = rescale(ref_fsimc, 0, 255);
    comp_fsimc = true;
end

[S, N] = size(Z);

measures.param2 = zeros(1, N);
measures.param1 = zeros(1, N);
measures.psnr = zeros(S, N);
measures.fsim = zeros(S,N);
measures.ssim = zeros(S,N);
measures.msam = zeros(S,N);
measures.rel_norm = zeros(S, N);

for k=1:S
    measures.param1(k) = Z{k, 1}.(param1_name);
    for j=1:N
        if ~isa(Z{k, j}, 'MException')
            measures.param2(j) = Z{k, j}.(param2_name);
        else
            measures.param2(j) = max(measures.param2(j), measures.param2(j-1) + eps);
        end

        for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
            n = i{1};
            measures.(n)(k, j) = get_measure(Z{k, j}, n);
        end
        
        if comp_fsimc
            if isa(Z{k, j}, 'MException')
                measures.fsimc(k, j) = 0;
            else
                [~, measures.fsimc(k, j)] = FeatureSIM(O, rescale(Z{k, j}.L, 0, 255));
            end
        end
        
    end
end

if comp_fsimc
    [best.fsimc.value, best.fsimc.index_param2] = max(measures.fsimc,[],2);
    [~, best.fsimc.index_param1] = max(best.fsimc.value);

    best.fsimc.index_param2 = best.fsimc.index_param2(best.fsimc.index_param1);
    best.fsimc.value = best.fsimc.value(best.fsimc.index_param1);
    
    best.fsimc.param1 = measures.param1(best.fsimc.index_param1);
    best.fsimc.param2 = measures.param2(best.fsimc.index_param2);
    best.fsimc.L = Z{best.fsimc.index_param1, best.fsimc.index_param2}.L;
    best.fsimc.E = Z{best.fsimc.index_param1, best.fsimc.index_param2}.E;
end
    
for i = {'psnr', 'fsim', 'ssim', 'rel_norm', 'msam'}
    n = i{1};
    if ~strcmp('msam', n) && ~strcmp('rel_norm', n)
        [best.(n).value, best.(n).index_param2] = max(measures.(n),[],2);
        [~, best.(n).index_param1] = max(best.(n).value);
    else
        [best.(n).value, best.(n).index_param2] = min(measures.(n),[],2);
        [~, best.(n).index_param1] = min(best.(n).value);
    end
    best.(n).index_param2 = best.(n).index_param2(best.(n).index_param1);
    best.(n).value = best.(n).value(best.(n).index_param1);
    
    best.(n).param1 = measures.param1(best.(n).index_param1);
    best.(n).param2 = measures.param2(best.(n).index_param2);
    best.(n).L = Z{best.(n).index_param1, best.(n).index_param2}.L;
    best.(n).E = Z{best.(n).index_param1, best.(n).index_param2}.E;
end

end