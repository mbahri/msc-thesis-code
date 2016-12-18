% [O, X] = yale_sp(1, 0.3);

k = 1;

f_psnr = figure;
[ha, pos] = tight_subplot(6, 3, [0.025 0.0], [0.0 0.025], [0.0 0.0]);
axes(ha(1)); ddisp(best.brtf.psnr.L(:,:,k)); title('BRTF');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(2)); ddisp(best.ten_rpca.psnr.L(:,:,k)); title('RPCA');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(3)); ddisp(best.ten_brpca.psnr.L(:,:,k)); title('BRPCA');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(4)); ddisp(best.ten_orpca.psnr.L(:,:,k)); title('RHOSVD');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(5)); ddisp(best.ten_rcpd.psnr.L(:,:,k)); title('RCPD');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(6)); ddisp(best.horpca_s.psnr.L(:,:,k)); title('HORPCA-S');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(7)); ddisp(best.horpca_s_tc.psnr.L(:,:,k)); title('HORPCA-S TC');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(8)); ddisp(best.rpca2d_l1.psnr.L(:,:,k)); title('RPCA2D L1');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(9)); ddisp(best.rpca2d_l2.psnr.L(:,:,k)); title('RPCA2D L2');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(10)); ddisp(best_2d.rpca2d_gl_l1.psnr.L(:,:,k)); title('RPCA2D GL L1');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(11)); ddisp(best_2d.rpca2d_gl_l2.psnr.L(:,:,k)); title('RPCA2D GL L2');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(12)); ddisp(best.cvpr2014_tsvd.psnr.L(:,:,k)); title('CVPR2014');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(13)); ddisp(best.cvpr2016_tnn.psnr.L(:,:,k)); title('CVPR2016');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(14)); ddisp(best_2d.cauchy_st.psnr.L(:,:,k)); title('Cauchy ST');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(15)); ddisp(best_2d.welsh_st.psnr.L(:,:,k)); title('Welsh ST');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(16)); ddisp(best_2d.nctrpca.psnr.L(:,:,k)); title('NC TRPCA');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(17)); ddisp(O(:,:,k)); title('Original');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(18)); ddisp(X(:,:,k)); title('Noisy');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);

f_psnr.Units = 'centimeters';
f_psnr.Position = [0 0 15 23];
f_psnr.PaperPositionMode = 'manual';
f_psnr.PaperUnits = 'centimeters';
f_psnr.PaperPosition = [0 0 15 23];
saveas(f_psnr, my_sprintf('grid_%s_psnr_%s_%0.2f.png', dataset, noise_type, noise_level), 'png');

k = 3;

f_psnr = figure;
[ha, pos] = tight_subplot(6, 3, [0.025 0.0], [0.0 0.025], [0.0 0.0]);
axes(ha(1)); ddisp(best.brtf.psnr.L(:,:,k)); title('BRTF'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(2)); ddisp(best.ten_rpca.psnr.L(:,:,k)); title('RPCA'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(3)); ddisp(best.ten_brpca.psnr.L(:,:,k)); title('BRPCA'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(4)); ddisp(best.ten_orpca.psnr.L(:,:,k)); title('RHOSVD'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(5)); ddisp(best.ten_rcpd.psnr.L(:,:,k)); title('RCPD'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(6)); ddisp(best.horpca_s.psnr.L(:,:,k)); title('HORPCA-S'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(7)); ddisp(best.horpca_s_tc.psnr.L(:,:,k)); title('HORPCA-S TC'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(8)); ddisp(best.rpca2d_l1.psnr.L(:,:,k)); title('RPCA2D L1'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(9)); ddisp(best.rpca2d_l2.psnr.L(:,:,k)); title('RPCA2D L2'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(10)); ddisp(best_2d.rpca2d_gl_l1.psnr.L(:,:,k)); title('RPCA2D GL L1'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(11)); ddisp(best_2d.rpca2d_gl_l2.psnr.L(:,:,k)); title('RPCA2D GL L2'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(12)); ddisp(best.cvpr2014_tsvd.psnr.L(:,:,k)); title('CVPR2014'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(13)); ddisp(best.cvpr2016_tnn.psnr.L(:,:,k)); title('CVPR2016'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(14)); ddisp(best_2d.cauchy_st.psnr.L(:,:,k)); title('Cauchy ST'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(15)); ddisp(best_2d.welsh_st.psnr.L(:,:,k)); title('Welsh ST'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(16)); ddisp(best_2d.nctrpca.psnr.L(:,:,k)); title('NC TRPCA'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(17)); ddisp(O(:,:,k)); title('Original'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(18)); ddisp(X(:,:,k)); title('Noisy'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);

f_psnr.Units = 'centimeters';
f_psnr.Position = [0 0 15 23];
f_psnr.PaperPositionMode = 'manual';
f_psnr.PaperUnits = 'centimeters';
f_psnr.PaperPosition = [0 0 15 23];
saveas(f_psnr, my_sprintf('grid_%s_psnr_%s_%0.2f_zoom_3.png', dataset, noise_type, noise_level), 'png');

k = 1;

f_fsim = figure;
[ha, pos] = tight_subplot(6, 3, [0.025 0.0], [0.0 0.025], [0.0 0.0]);
axes(ha(1)); ddisp(best.brtf.fsim.L(:,:,k)); title('BRTF');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(2)); ddisp(best.ten_rpca.fsim.L(:,:,k)); title('RPCA');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(3)); ddisp(best.ten_brpca.fsim.L(:,:,k)); title('BRPCA');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(4)); ddisp(best.ten_orpca.fsim.L(:,:,k)); title('RHOSVD');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(5)); ddisp(best.ten_rcpd.fsim.L(:,:,k)); title('RCPD'); %xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(6)); ddisp(best.horpca_s.fsim.L(:,:,k)); title('HORPCA-S');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(7)); ddisp(best.horpca_s_tc.fsim.L(:,:,k)); title('HORPCA-S TC'); %xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(8)); ddisp(best.rpca2d_l1.fsim.L(:,:,k)); title('RPCA2D L1'); %xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(9)); ddisp(best.rpca2d_l2.fsim.L(:,:,k)); title('RPCA2D L2'); %xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(10)); ddisp(best_2d.rpca2d_gl_l1.fsim.L(:,:,k)); title('RPCA2D GL L1');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(11)); ddisp(best_2d.rpca2d_gl_l2.fsim.L(:,:,k)); title('RPCA2D GL L2');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(12)); ddisp(best.cvpr2014_tsvd.fsim.L(:,:,k)); title('CVPR2014');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(13)); ddisp(best.cvpr2016_tnn.fsim.L(:,:,k)); title('CVPR2016');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(14)); ddisp(best_2d.cauchy_st.fsim.L(:,:,k)); title('Cauchy ST'); %xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(15)); ddisp(best_2d.welsh_st.fsim.L(:,:,k)); title('Welsh ST');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(16)); ddisp(best_2d.nctrpca.fsim.L(:,:,k)); title('NC TRPCA');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(17)); ddisp(O(:,:,k)); title('Original'); %xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(18)); ddisp(X(:,:,k)); title('Noisy');% xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);

f_fsim.Units = 'centimeters';
f_fsim.Position = [0 0 15 23];
f_fsim.PaperPositionMode = 'manual';
f_fsim.PaperUnits = 'centimeters';
f_fsim.PaperPosition = [0 0 15 23];
saveas(f_fsim, my_sprintf('grid_%s_fsim_%s_%0.2f.png', dataset, noise_type, noise_level), 'png');

k = 3;

f_fsim = figure;
[ha, pos] = tight_subplot(6, 3, [0.025 0.0], [0.0 0.025], [0.0 0.0]);
axes(ha(1)); ddisp(best.brtf.fsim.L(:,:,k)); title('BRTF'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(2)); ddisp(best.ten_rpca.fsim.L(:,:,k)); title('RPCA'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(3)); ddisp(best.ten_brpca.fsim.L(:,:,k)); title('BRPCA'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(4)); ddisp(best.ten_orpca.fsim.L(:,:,k)); title('RHOSVD'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(5)); ddisp(best.ten_rcpd.fsim.L(:,:,k)); title('RCPD'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(6)); ddisp(best.horpca_s.fsim.L(:,:,k)); title('HORPCA-S'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(7)); ddisp(best.horpca_s_tc.fsim.L(:,:,k)); title('HORPCA-S TC'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(8)); ddisp(best.rpca2d_l1.fsim.L(:,:,k)); title('RPCA2D L1'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(9)); ddisp(best.rpca2d_l2.fsim.L(:,:,k)); title('RPCA2D L2'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(10)); ddisp(best_2d.rpca2d_gl_l1.fsim.L(:,:,k)); title('RPCA2D GL L1'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(11)); ddisp(best_2d.rpca2d_gl_l2.fsim.L(:,:,k)); title('RPCA2D GL L2'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(12)); ddisp(best.cvpr2014_tsvd.fsim.L(:,:,k)); title('CVPR2014'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(13)); ddisp(best.cvpr2016_tnn.fsim.L(:,:,k)); title('CVPR2016'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(14)); ddisp(best_2d.cauchy_st.fsim.L(:,:,k)); title('Cauchy ST'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(15)); ddisp(best_2d.welsh_st.fsim.L(:,:,k)); title('Welsh ST'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(16)); ddisp(best_2d.nctrpca.fsim.L(:,:,k)); title('NC TRPCA'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(17)); ddisp(O(:,:,k)); title('Original'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);
axes(ha(18)); ddisp(X(:,:,k)); title('Noisy'); xlim([0.5000   84.5000]); ylim([23.9435  119.9435]);

f_fsim.Units = 'centimeters';
f_fsim.Position = [0 0 15 23];
f_fsim.PaperPositionMode = 'manual';
f_fsim.PaperUnits = 'centimeters';
f_fsim.PaperPosition = [0 0 15 23];
saveas(f_fsim, my_sprintf('grid_%s_fsim_%s_%0.2f_zoom_3.png', dataset, noise_type, noise_level), 'png');