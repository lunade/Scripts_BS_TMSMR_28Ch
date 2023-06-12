% Script to show SNR comp
% 
% Lucia Navarro de Lara, PhD
% MGH, Martinos Center
% 04/25/2023
%
% Instructions:
% PROCESSING FOLDER SHOULD HAVE PROCESSED DATA from calc_SNR_g_factors_invivo_TMSMR_20Ch/calc_SNR_g_factors_invivo_TMSMR_32Ch
% mask_manual2.mat is used to mask the maps
% FIGURE_FOLDER is the folder to save figures
%clear all;

%figure_folder='/autofs/space/tmsbox_002/users/lucia/figures/paper_BS_RFCoil/revision'
figure_folder='FIGURE_FOLDER'
base_proce_folder='PROCSSEING_FOLDER/'
load([base_proce_folder,'/snr_mapping_TMSMR_32Ch_comp_July_2022_nonorm.mat'])
load([base_proce_folder,'/mask_manual2.mat'])
sli_to_plot = 20:2:58;
sli_to_average=20:58;
dim = numel(snr_with_tms_norm(:,1,1));
x = linspace(-dim/2,dim/2,dim)
mask=msk_temp(:,:,5:end-4);
%%
% mask brain to get rid of BET;
snr_32Ch_masked=snr_32Ch_final.*mask;
snr_with_tms_only_brain=snr_with_tms_norm.*mask;


crop_x=20;
crop_y=20;
% %
 h_1=figure(28),imagesc(vol2mos(snr_32Ch_masked(crop_x:end-crop_x,crop_y:end-crop_y,sli_to_plot)),[0 200]),colormap(jet),axis image off;
 saveas(h_1,[figure_folder,'/SNR_32Ch_v7_nonorm.tif']);
 h_2=figure(29),imagesc(vol2mos(snr_with_tms_only_brain(crop_x:end-crop_x,crop_y:end-crop_y,sli_to_plot)),[0 200]),colormap(jet),axis image off;
 saveas(h_2,[figure_folder,'/SNR_TMSMR_v7_nonorm.tif']);
 for i=1:numel(sli_to_average)
     snr_gain_32Ch_slices_mean(i)=mean(snr_gain_32Ch(:,:,sli_to_average(i)),'all','omitNaN');
 end
 save([base_proce_folder,'/snr_gain_TMSvs32_v7_no_norm.mat'], "snr_gain_32Ch_slices_mean");
 
