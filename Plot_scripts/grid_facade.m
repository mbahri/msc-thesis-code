% [O, X] = facade_sp(1, 0.3);

figure
subplot(6,3,1); ddisp(best.brtf.psnr.L); title('BRTF');
subplot(6,3,2); ddisp(best.ten_rpca.psnr.L); title('RPCA');
subplot(6,3,3); ddisp(best.ten_brpca.psnr.L); title('BRPCA');
subplot(6,3,4); ddisp(best.ten_orpca.psnr.L); title('HORSVD');
subplot(6,3,5); ddisp(best.ten_rcpd.psnr.L); title('RCPD');
subplot(6,3,6); ddisp(best.horpca_s.psnr.L); title('HORPCA-S');
subplot(6,3,7); ddisp(best.horpca_s_tc.psnr.L); title('HORPCA-S TC');
subplot(6,3,8); ddisp(best.rpca2d_l1.psnr.L); title('RPCA2D L1');
subplot(6,3,9); ddisp(best.rpca2d_l2.psnr.L); title('RPCA2D L2');
subplot(6,3,10); ddisp(best_2d.rpca2d_gl_l1.psnr.L); title('RPCA2D GL L1');
subplot(6,3,11); ddisp(best_2d.rpca2d_gl_l2.psnr.L); title('RPCA2D GL L2');
subplot(6,3,12); ddisp(best.cvpr2014_tsvd.psnr.L); title('CVPR2014');
subplot(6,3,13); ddisp(best.cvpr2016_tnn.psnr.L); title('CVPR2016');
subplot(6,3,14); ddisp(best_2d.cauchy_st.psnr.L); title('Cauchy ST');
subplot(6,3,15); ddisp(best_2d.welsh_st.psnr.L); title('Welsh ST');
subplot(6,3,16); ddisp(best_2d.nctrpca.psnr.L); title('NCTRPCA');
subplot(6,3,17); ddisp(O); title('Original');
subplot(6,3,18); ddisp(X); title('Noisy');

figure
subplot(6,3,1); ddisp(best.brtf.fsimc.L); title('BRTF');
subplot(6,3,2); ddisp(best.ten_rpca.fsimc.L); title('RPCA');
subplot(6,3,3); ddisp(best.ten_brpca.fsimc.L); title('BRPCA');
subplot(6,3,4); ddisp(best.ten_orpca.fsimc.L); title('HORSVD');
subplot(6,3,5); ddisp(best.ten_rcpd.fsimc.L); title('RCPD');
subplot(6,3,6); ddisp(best.horpca_s.fsimc.L); title('HORPCA-S');
subplot(6,3,7); ddisp(best.horpca_s_tc.fsimc.L); title('HORPCA-S TC');
subplot(6,3,8); ddisp(best.rpca2d_l1.fsimc.L); title('RPCA2D L1');
subplot(6,3,9); ddisp(best.rpca2d_l2.fsimc.L); title('RPCA2D L2');
subplot(6,3,10); ddisp(best_2d.rpca2d_gl_l1.fsimc.L); title('RPCA2D GL L1');
subplot(6,3,11); ddisp(best_2d.rpca2d_gl_l2.fsimc.L); title('RPCA2D GL L2');
subplot(6,3,12); ddisp(best.cvpr2014_tsvd.fsimc.L); title('CVPR2014');
subplot(6,3,13); ddisp(best.cvpr2016_tnn.fsimc.L); title('CVPR2016');
subplot(6,3,14); ddisp(best_2d.cauchy_st.fsimc.L); title('Cauchy ST');
subplot(6,3,15); ddisp(best_2d.welsh_st.fsimc.L); title('Welsh ST');
subplot(6,3,16); ddisp(best_2d.nctrpca.fsimc.L); title('NCTRPCA');
subplot(6,3,17); ddisp(O); title('Original');
subplot(6,3,18); ddisp(X); title('Noisy');

figure
subplot(6,3,1); ddisp(best.brtf.msam.L); title('BRTF');
subplot(6,3,2); ddisp(best.ten_rpca.msam.L); title('RPCA');
subplot(6,3,3); ddisp(best.ten_brpca.msam.L); title('BRPCA');
subplot(6,3,4); ddisp(best.ten_orpca.msam.L); title('HORSVD');
subplot(6,3,5); ddisp(best.ten_rcpd.msam.L); title('RCPD');
subplot(6,3,6); ddisp(best.horpca_s.msam.L); title('HORPCA-S');
subplot(6,3,7); ddisp(best.horpca_s_tc.msam.L); title('HORPCA-S TC');
subplot(6,3,8); ddisp(best.rpca2d_l1.msam.L); title('RPCA2D L1');
subplot(6,3,9); ddisp(best.rpca2d_l2.msam.L); title('RPCA2D L2');
subplot(6,3,10); ddisp(best_2d.rpca2d_gl_l1.msam.L); title('RPCA2D GL L1');
subplot(6,3,11); ddisp(best_2d.rpca2d_gl_l2.msam.L); title('RPCA2D GL L2');
subplot(6,3,12); ddisp(best.cvpr2014_tsvd.msam.L); title('CVPR2014');
subplot(6,3,13); ddisp(best.cvpr2016_tnn.msam.L); title('CVPR2016');
subplot(6,3,14); ddisp(best_2d.cauchy_st.msam.L); title('Cauchy ST');
subplot(6,3,15); ddisp(best_2d.welsh_st.msam.L); title('Welsh ST');
subplot(6,3,16); ddisp(best_2d.nctrpca.msam.L); title('NCTRPCA');
subplot(6,3,17); ddisp(O); title('Original');
subplot(6,3,18); ddisp(X); title('Noisy');