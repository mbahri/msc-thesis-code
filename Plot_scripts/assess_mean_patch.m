clear;

subfolder = 'Yale_Patch';
dataset = 'yale';
noise_type = 'patch';
noise_level = 160;

[O, ~] = yale_patch(1, noise_level);
O = O(:,:,1:64);
True_M = mean(O, 3);

% fprintf('Loading Mehdi...\n');
% method = 'mehdi';
% to_load = sprintf('/vol/bitbucket/mb2215/Thesis/Results/%s/results_%s_%s_%s_%f.mat', subfolder, dataset, method, noise_type, noise_level);
% load(to_load)

% LOL

load('/vol/bitbucket/mb2215/Thesis/Results/Yale_Patch/results_yale_mehdi_gl_l1_patch_160.000000.mat')
[N1, N2] = size(results_mehdi.rpca2d_gl_l1);
for i = 1:N1
    for j = 1:N2
        scores.rpca2d_gl_l1{i, j} = MSIQA(True_M, results_mehdi.rpca2d_gl_l1{i, j}.M);
        psnr.rpca2d_gl_l1(i, j) = scores.rpca2d_gl_l1{i, j}.psnr;
        fsim.rpca2d_gl_l1(i, j) = scores.rpca2d_gl_l1{i, j}.fsim;
        ssim.rpca2d_gl_l1(i, j) = scores.rpca2d_gl_l1{i, j}.ssim;
        msam.rpca2d_gl_l1(i, j) = scores.rpca2d_gl_l1{i, j}.msam;
        rel_norm.rpca2d_gl_l1(i, j) = scores.rpca2d_gl_l1{i, j}.rel_norm;
    end
end

% [mpsnr, impsnr] = max(psnr.rpca2d_gl_l1);
% [mmpsnr, jmpsnr] = max(mpsnr);
% toplot = psnr.rpca2d_gl_l1(:,jmpsnr);
% impsnr = impsnr(jmpsnr);
% 
% plot(toplot);
% 
% ddisp(results_mehdi.rpca2d_gl_l1{impsnr, jmpsnr}.M)

[mfsim, imfsim] = max(fsim.rpca2d_gl_l1);
[mmfsim, jmfsim] = max(mfsim);
toplot = fsim.rpca2d_gl_l1(:,jmfsim);
imfsim = imfsim(jmfsim);

plot(toplot)

ddisp(results_mehdi.rpca2d_gl_l1{imfsim, jmfsim}.M)

scores.rpca2d_gl_l1{imfsim, jmfsim}

% [mrel_norm, imrel_norm] = max(rel_norm.rpca2d_gl_l1);
% [mmrel_norm, jmrel_norm] = max(mrel_norm);
% toplot = rel_norm.rpca2d_gl_l1(:,jmrel_norm);
% imrel_norm = imrel_norm(jmrel_norm);
% 
% plot(toplot)
% 
% ddisp(results_mehdi.rpca2d_gl_l1{imrel_norm, jmrel_norm}.M)

load('/vol/bitbucket/mb2215/Thesis/Results/Yale_Patch/results_yale_mehdi_gl_l2_patch_160.000000.mat')
[N1, N2] = size(results_mehdi.rpca2d_gl_l2);
for i = 1:N1
    for j = 1:N2
        scores.rpca2d_gl_l2{i, j} = MSIQA(True_M, results_mehdi.rpca2d_gl_l2{i, j}.M);
        psnr.rpca2d_gl_l2(i, j) = scores.rpca2d_gl_l2{i, j}.psnr;
        fsim.rpca2d_gl_l2(i, j) = scores.rpca2d_gl_l2{i, j}.fsim;
        ssim.rpca2d_gl_l2(i, j) = scores.rpca2d_gl_l2{i, j}.ssim;
        msam.rpca2d_gl_l2(i, j) = scores.rpca2d_gl_l2{i, j}.msam;
        rel_norm.rpca2d_gl_l2(i, j) = scores.rpca2d_gl_l2{i, j}.rel_norm;
    end
end

% [mpsnr, impsnr] = max(psnr.rpca2d_gl_l2);
% [mmpsnr, jmpsnr] = max(mpsnr);
% toplot = psnr.rpca2d_gl_l2(:,jmpsnr);
% impsnr = impsnr(jmpsnr);
% 
% plot(toplot);
% 
% ddisp(results_mehdi.rpca2d_gl_l2{impsnr, jmpsnr}.M)

[mfsim, imfsim] = max(fsim.rpca2d_gl_l2);
[mmfsim, jmfsim] = max(mfsim);
toplot = fsim.rpca2d_gl_l2(:,jmfsim);
imfsim = imfsim(jmfsim);

plot(toplot)

ddisp(results_mehdi.rpca2d_gl_l2{imfsim, jmfsim}.M)

scores.rpca2d_gl_l2{imfsim, jmfsim}

% [mrel_norm, imrel_norm] = max(rel_norm.rpca2d_gl_l2);
% [mmrel_norm, jmrel_norm] = max(mrel_norm);
% toplot = rel_norm.rpca2d_gl_l2(:,jmrel_norm);
% imrel_norm = imrel_norm(jmrel_norm);
% 
% plot(toplot)
% 
% ddisp(results_mehdi.rpca2d_gl_l2{imrel_norm, jmrel_norm}.M)