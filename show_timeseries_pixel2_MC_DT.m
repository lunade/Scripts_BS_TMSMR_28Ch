% Script to show SNR comp
% 
% Lucia Navarro de Lara, PhD
% MGH, Martinos Center
% 04/25/2023
%
% Instructions:
% PROCESSING FOLDER SHOULD HAVE FMRI PROCESSED DATA 
% mask_manual2.mat is used to mask the maps
% FIGURE_FOLDER is the folder to save figures


epi_with_tms_filepath = 'PROCESSING_FOLDER';

nii_epi_tmsmr28ch_mc_dt= load_nifti([epi_with_tms_filepath,'epi_tmsmr28ch_mc_dt.nii.gz']);
nii_epi_mc_dt=rot90(nii_epi_tmsmr28ch_mc_dt.vol);

nii_epi_tmsmr28ch= load_nifti([epi_with_tms_filepath,'epi_tmsmr28ch.nii.gz']);
nii_epi=rot90(nii_epi_tmsmr28ch.vol);
% voxel 1 slice 22, index(35 65)
%voxel 2 slice 24, index(23,55)
%voxel 3 slice 26, index(36,37)
% Plot temporal signal
f_b=figure(1);
colormap(gray(256))
c_min=70; c_max=830;
imagesc(squeeze(nii_epi_mc_dt(:,:,24,21)));
daspect([1 1 1]);axis off;caxis([c_min c_max]);
f_sig=figure(2);
signal=squeeze(nii_epi_mc_dt(23,55,24,21:end));

plot(signal,'r','LineWidth',4); hold on
axis_c=gca;
axis_c.GridAlphaMode='manual';
axis_c.GridAlpha=0.5;
axis_c.GridLineStyle='--';
axis_c.LineWidth=2;
axis_c.GridColor='k';
axis_c.YGrid='on';
axis_c.XGrid='off';
axis_c.FontSize=20;



signal_mean=mean(signal);
noise=std(signal);
tSNR=signal_mean/noise;
axis([1 180 signal_mean-(0.25*signal_mean)  signal_mean+(0.25*signal_mean)]);
%%
fig_folder='/autofs/space/tmsbox_002/users/lucia/figures/ISMRM_2022/';
fig_folder='FIGURES_FOLDER';
%%
saveas(f_sig,[fig_folder,'EPI_brain_1_timesignal_mc_dt.tif']);
