
% [O, X] = yale_patch(1, 160);
% [O, X] = yale_sp(1, 0.3);

figure
k = 1;
subplot(6,3,1); ddisp(best.brtf.psnr.L(:,:,k)); title('BRTF');
subplot(6,3,2); ddisp(best.ten_rpca.psnr.L(:,:,k)); title('RPCA');
subplot(6,3,3); ddisp(best.ten_brpca.psnr.L(:,:,k)); title('BRPCA');
subplot(6,3,4); ddisp(best.ten_orpca.psnr.L(:,:,k)); title('HORSVD');
subplot(6,3,5); ddisp(best.ten_rcpd.psnr.L(:,:,k)); title('RCPD');
subplot(6,3,6); ddisp(best.horpca_s.psnr.L(:,:,k)); title('HORPCA-S');
subplot(6,3,7); ddisp(best.horpca_s_tc.psnr.L(:,:,k)); title('HORPCA-S TC');
subplot(6,3,8); ddisp(best.rpca2d_l1.psnr.L(:,:,k)); title('RPCA2D L1');
subplot(6,3,9); ddisp(best.rpca2d_l2.psnr.L(:,:,k)); title('RPCA2D L2');
subplot(6,3,10); ddisp(best_2d.rpca2d_gl_l1.psnr.L(:,:,k)); title('RPCA2D GL L1');
subplot(6,3,11); ddisp(best_2d.rpca2d_gl_l2.psnr.L(:,:,k)); title('RPCA2D GL L2');
subplot(6,3,12); ddisp(best.cvpr2014_tsvd.psnr.L(:,:,k)); title('CVPR2014');
subplot(6,3,13); ddisp(best.cvpr2016_tnn.psnr.L(:,:,k)); title('CVPR2016');
subplot(6,3,14); ddisp(best_2d.cauchy_st.psnr.L(:,:,k)); title('Cauchy ST');
subplot(6,3,15); ddisp(best_2d.welsh_st.psnr.L(:,:,k)); title('Welsh ST');
subplot(6,3,16); ddisp(best_2d.nctrpca.psnr.L(:,:,k)); title('NCTRPCA');
subplot(6,3,17); ddisp(O(:,:,k)); title('Original');
subplot(6,3,18); ddisp(X(:,:,k)); title('Noisy');

figure
subplot(6,3,1); ddisp(best.brtf.fsim.L(:,:,k)); title('BRTF');
subplot(6,3,2); ddisp(best.ten_rpca.fsim.L(:,:,k)); title('RPCA');
subplot(6,3,3); ddisp(best.ten_brpca.fsim.L(:,:,k)); title('BRPCA');
subplot(6,3,4); ddisp(best.ten_orpca.fsim.L(:,:,k)); title('HORSVD');
subplot(6,3,5); ddisp(best.ten_rcpd.fsim.L(:,:,k)); title('RCPD');
subplot(6,3,6); ddisp(best.horpca_s.fsim.L(:,:,k)); title('HORPCA-S');
subplot(6,3,7); ddisp(best.horpca_s_tc.fsim.L(:,:,k)); title('HORPCA-S TC');
subplot(6,3,8); ddisp(best.rpca2d_l1.fsim.L(:,:,k)); title('RPCA2D L1');
subplot(6,3,9); ddisp(best.rpca2d_l2.fsim.L(:,:,k)); title('RPCA2D L2');
subplot(6,3,10); ddisp(best_2d.rpca2d_gl_l1.fsim.L(:,:,k)); title('RPCA2D GL L1');
subplot(6,3,11); ddisp(best_2d.rpca2d_gl_l2.fsim.L(:,:,k)); title('RPCA2D GL L2');
subplot(6,3,12); ddisp(best.cvpr2014_tsvd.fsim.L(:,:,k)); title('CVPR2014');
subplot(6,3,13); ddisp(best.cvpr2016_tnn.fsim.L(:,:,k)); title('CVPR2016');
subplot(6,3,14); ddisp(best_2d.cauchy_st.fsim.L(:,:,k)); title('Cauchy ST');
subplot(6,3,15); ddisp(best_2d.welsh_st.fsim.L(:,:,k)); title('Welsh ST');
subplot(6,3,16); ddisp(best_2d.nctrpca.fsim.L(:,:,k)); title('NCTRPCA');
subplot(6,3,17); ddisp(O(:,:,k)); title('Original');
subplot(6,3,18); ddisp(X(:,:,k)); title('Noisy');