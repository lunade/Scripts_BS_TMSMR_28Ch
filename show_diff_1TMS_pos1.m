
% Show difference in SNR maps 
% Lucia Navarro de Lara, PhD
% MGH, Martinos Center
% 04/25/2023
%
% Instructions:
% PROCESSING FOLDER SHOULD HAVE PROCSSED DATA from SNR_norm_TMSMR_1TMS_pos1
% colormap_dif.mat to use same colormap
% FIGURE_FOLDER is the folder to save figures




%%load data
pathrootprocess = 'PROCESSING FOLDER';
load([pathrootprocess, '/results_norm_pos1.mat']);
load([pathrootprocess, '/colormap_dif.mat']);
figure_folder='FIGURE FOLDER'
%Plot coronal normalized data
%%Display diff maps snr_cov
sli_to_plot = 11:1:22;
sli_to_plot_2= 20:2:80;
sli_to_plot_3= 20:2:96;
slc_sag=66;slc_tra=24;slc_cor=34;
colormap(cmap_france);
snr_cov_norm_masked=snr_cov_dif_norm.*msk;
snr_cov_norm_masked_perm=permute(snr_cov_norm_masked,[3,1,2]);
snr_cov_norm_masked_perm2=permute(snr_cov_norm_masked,[3,2,1]);
snr_cov_no_tms_norm_masked=snr_cov_no_tms_norm.*msk;
snr_cov_no_tms_norm_masked_perm=permute(snr_cov_no_tms_norm_masked,[3,1,2]);
snr_cov_no_tms_norm_masked_perm2=permute(snr_cov_no_tms_norm_masked,[3,2,1]);

h=figure(6),imagesc([(vol2mos(rot90(snr_cov_norm_masked(:,:,sli_to_plot),1)))],[-40 40]),axis image off, colormap(cmap_france);
saveas(h,[figure_folder,'/snr_diff_measured_pos1_Nov.tif'])

h=figure(7),imagesc([(vol2mos(rot90(snr_cov_no_tms_norm_masked(:,:,sli_to_plot),1)))],[0 1200]),axis image off, colormap("jet");
saveas(h,[figure_folder,'/snr_no_tms_measured_pos1_Nov.tif'])

