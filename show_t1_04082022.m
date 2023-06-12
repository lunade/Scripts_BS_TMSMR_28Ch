
%Script to show T1/MPRAGE
% 
% Lucia Navarro de Lara, PhD
% MGH, Martinos Center
% 04/25/2023
%
% Instructions:
% BASE FOLDER SHOULD HAVE PROCESSED DATA from calc_SNR_g_factors_invivo_TMSMR_20Ch/calc_SNR_g_factors_invivo_TMSMR_32Ch
% mask_manual2.mat is used to mask the maps
% FIGURE_FOLDER is the folder to save figures


base_folder = '/autofs/cluster/tmslab/lucia/TMSMR_28ch_SNR/04082022/DICOM/T1_MEMPRAGE_SAG_NONORM_RMS_0005';
%t1_with_tms_filepath = [base_folder,'/DICOM/T1_MEMPRAGE_SAG_NONORM_RMS_0005'];
 
%threshold=3000;
nii_t1_tmsmr28ch= load_nifti([base_folder,'/t1_tmsmr28Ch.nii']);
nii_t1=rot90(nii_t1_tmsmr28ch.vol,1);
%nii_t1_ori=flipud(nii_t1);
%%

c_min=0; c_max=264;
%Sag
F_sag=figure;
imagesc(squeeze(nii_t1(:,:,88)));
daspect([1 1 1]);
axis off;caxis([c_min c_max]);colormap(gray(256))
% Cor flips to have radiological ori
F_cor=figure;
imagesc(squeeze(flip(nii_t1(:,128,:),3)));
daspect([1 1 1]);axis off;caxis([c_min c_max]);colormap(gray(256))
% Tra flips to have radiological ori
F_tra=figure;
imagesc(squeeze(flip(nii_t1(80,:,:),3)));
daspect([1 1 1]);axis off;caxis([c_min c_max]);colormap(gray(256))

%%
cd '/autofs/space/tmsbox_002/users/lucia/figures/TMSMR28CH_Paper';
saveas(F_sag,'T1_TMSMRMR28Ch_sag_v2.tif');
saveas(F_cor,'T1_TMSMRMR28Ch_cor_v2.tif');
saveas(F_tra,'T1_TMSMRMR28Ch_tra_v2.tif');