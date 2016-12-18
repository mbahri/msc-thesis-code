% load_highway;
% GT_frames = [17 105];
% O = X(:,:,GT_frames);
% GT = GT(:,:,GT_frames);

load_hall;
O = O(:,:,GT_frames);

dataset = 'hall';

best_r_idx('rpca') = 1;

algorithms = {'rpca',  'brpca', 'orpca', 'rcpd_sub', 'rcpd_lin' ...
              'rpca2d_l1', 'rpca2d_l2', 'rpca2d_gl_l1', 'rpca2d_gl_l2', ...
              'tnn', 'brtf', 'horpca_s', ...
              'tsvd', 'welsh_st', 'cauchy_st', 'nctrpca'};
          
names = {'Ten RPCA', 'Ten BRPCA', 'RHOVSD', 'RCPD Sub', 'RCPD Lin', ...
         '2D L1', '2D L2', '2D GL L1', '2D GL L2', ...
         'CVPR2016', 'BRTF', 'HORPCA-S', 'CVPR2014', 'Welsh ST', ...
         'Cauchy ST', 'NCTRPCA'};
          
n_alg = length(algorithms);

part1 = figure;
[ha, ~] = tight_subplot(1 + n_alg/2, 4, [0.0 0.01], [0.01 0.01], [0.05 0.0]);
j = 0;

% Header with references
axes(ha(1)); ddisp(O(:,:,1)); yl = ylabel('Reference'); yl.Units = 'Normalized'; yl.Position = [-0.05, 0.5, 0];
axes(ha(3)); ddisp(O(:,:,2)); 
axes(ha(2)); ddisp(GT(:,:,1));
axes(ha(4)); ddisp(GT(:,:,2));

for i=5:4:2*n_alg + 4
    j = j + 1;
    n = algorithms{j};
    
    if strcmp('rcpd_sub', n)
        bg = double(results.rcpd.sub.A{best_l_idx(n), best_r_idx(n)});
    elseif strcmp('rcpd_lin', n)
        bg = double(results.rcpd.lin.A{best_l_idx(n), best_r_idx(n)});
    else
        bg = double(results.(n).A{best_l_idx(n), best_r_idx(n)});
    end
%     bg = bg(:,:,GT_frames);
    
    fg = double(best_foreground(n));
%     fg = fg(:,:,GT_frames);
%     axes(ha(i)); ddisp(...
%         [...
%         bg(:,:,1) ...
%         bg(:,:,2) ...
%         fg(:,:,1) ...
%         fg(:,:,2) ...
%         ]...
%         ); ylabel(names{i});
    axes(ha(i)); ddisp(bg(:,:,1)); yl = ylabel(names{j}); yl.Units = 'Normalized'; yl.Position = [-0.05, 0.5, 0];
    axes(ha(i+2)); ddisp(bg(:,:,2)); 
    axes(ha(i+1)); ddisp(fg(:,:,1));
    axes(ha(i+3)); ddisp(fg(:,:,2));
end

part1.Units = 'centimeters';
part1.Position = [0 0 12.57 23];
% part1.Position = [0 0 14.34 23];
part1.PaperPositionMode = 'manual';
part1.PaperUnits = 'centimeters';
part1.PaperPosition = [0 0 12.57 23];
% part1.PaperPosition = [0 0 14.34 23];
saveas(part1, my_sprintf('part1_%s.png', dataset), 'png');

% figure
% [ha, ~] = tight_subplot(n_alg/2, 1, [0.01 0.0], [0.025 0.025], [0.0 0.0]);
% for j=n_alg/2 + 1:n_alg
%     i = j - (n_alg/2);
%     n = algorithms{i};
%     
%     if strcmp('rcpd_sub', n)
%         bg = double(results.rcpd.sub.A{best_l_idx(n), best_r_idx(n)});
%     elseif strcmp('rcpd_lin', n)
%         bg = double(results.rcpd.lin.A{best_l_idx(n), best_r_idx(n)});
%     else
%         bg = double(results.(n).A{best_l_idx(n), best_r_idx(n)});
%     end
%     
%     fg = double(best_foreground(n));
%     axes(ha(i)); ddisp(...
%         [...
%         bg(:,:,1) ...
%         bg(:,:,2) ...
%         fg(:,:,1) ...
%         fg(:,:,2) ...
%         ]...
%         ); ylabel(names{i});
% end

part2 = figure;
[ha, ~] = tight_subplot(1 + n_alg/2, 4, [0.0 0.01], [0.01 0.01], [0.05 0.0]);
% j = 0;

% Header with references
axes(ha(1)); ddisp(O(:,:,1)); yl = ylabel('Reference'); yl.Units = 'Normalized'; yl.Position = [-0.05, 0.5, 0];
axes(ha(3)); ddisp(O(:,:,2)); 
axes(ha(2)); ddisp(GT(:,:,1));
axes(ha(4)); ddisp(GT(:,:,2));

for k=2*n_alg + 9:4:4*n_alg + 8
    i = k - (2*n_alg + 4);
    j = j + 1;
    n = algorithms{j};
    
    if strcmp('rcpd_sub', n)
        bg = double(results.rcpd.sub.A{best_l_idx(n), best_r_idx(n)});
    elseif strcmp('rcpd_lin', n)
        bg = double(results.rcpd.lin.A{best_l_idx(n), best_r_idx(n)});
    else
        bg = double(results.(n).A{best_l_idx(n), best_r_idx(n)});
    end
%     bg = bg(:,:,GT_frames);
    
    fg = double(best_foreground(n));
%     fg = fg(:,:,GT_frames);
%     axes(ha(i)); ddisp(...
%         [...
%         bg(:,:,1) ...
%         bg(:,:,2) ...
%         fg(:,:,1) ...
%         fg(:,:,2) ...
%         ]...
%         ); ylabel(names{i});
    axes(ha(i)); ddisp(bg(:,:,1)); yl = ylabel(names{j}); yl.Units = 'Normalized'; yl.Position = [-0.05, 0.5, 0];
    axes(ha(i+2)); ddisp(bg(:,:,2)); 
    axes(ha(i+1)); ddisp(fg(:,:,1));
    axes(ha(i+3)); ddisp(fg(:,:,2));
end

part2.Units = 'centimeters';
part2.Position = [0 0 12.57 23];
% part2.Position = [0 0 14.34 23];
part2.PaperPositionMode = 'manual';
part2.PaperUnits = 'centimeters';
part2.PaperPosition = [0 0 12.57 23];
% part2.PaperPosition = [0 0 14.34 23];

saveas(part2, my_sprintf('part2_%s.png', dataset), 'png');