% load_highway;
% GT_frames = [17 105];
% O = X(:,:,GT_frames);
% GT = GT(:,:,GT_frames);

load_hall;
O = O(:,:,GT_frames);

dataset = 'hall';

algorithms = {'rpca2d_l1', 'tnn', 'horpca_s', ...
              'tsvd', 'welsh_st', 'cauchy_st', 'nctrpca'};
          
names = {'2D L1', 'CVPR2016', 'HORPCA-S', 'CVPR2014', 'Welsh ST', ...
         'Cauchy ST', 'NCTRPCA'};
          
n_alg = length(algorithms);

for j=1:n_alg
    n = algorithms{j};
    bg = double(results.(n).A{best_l_idx(n), best_r_idx(n)});
    fg = double(best_foreground(n));
    imwrite(bg(:,:,1), sprintf('%s_bg_1.png', n));
    imwrite(bg(:,:,2), sprintf('%s_bg_2.png', n));
    imwrite(fg(:,:,1), sprintf('%s_fg_1.png', n));
    imwrite(fg(:,:,2), sprintf('%s_fg_2.png', n));
end

