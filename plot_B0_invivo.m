
%Calculate and plot B0 maps
%Lucia Navarro de Lara, PhD
%MGH, Martinos Center
% 04/25/2023 lnavarrodelara@mgh.harvard.edu
%DICOM PATH in maps_dir
%delta_TE is information from the sequence
%
clear all, close all
curr_path = fileparts(mfilename('fullpath'))
cd(curr_path);
addpath(curr_path);
% Path to DICOM images of B0 map
maps_dir= B0_MAPS_DICOM;
%maps_dir='/cluster/tmslab/lucia/TMSMR_28ch_SNR/04082022/DICOM/GRE_FIELD_MAPPING_B0_0012'
% Path to a mask of the images
load('/autofs/space/tmsbox_002/users/lucia/Processing_TMSMR_28Ch/lost&found/mask_manual_b0.mat');
load('/autofs/space/tmsbox_002/users/lucia/Processing_TMSMR_28Ch/B0_maps/fine_b0_msk.mat');
% Parameters of the sequence
delta_TE=2.46e-3;


nii_b0=dicom_fm_import(maps_dir);
% Conversion to Hz
nii_b0_maps= pi*((double(nii_b0)-2048)/2048)/delta_TE/2/pi;
% To create contrast
msk_con=-200*(~msk);

%% Plot results
cd /autofs/space/tmsbox_002/users/lucia/figures/TMSMR28CH_Paper/
crop=18;
nii_b0_2display=(nii_b0_maps.*msk)+msk_con;
    fig1=figure;plot_range = [-100 100]; plot_range_dif=[-10 10];
    imagesc([vol2mos(nii_b0_2display(crop:end-crop,crop:end-crop,34:58))],plot_range),axis image, axis off,colorbar
    colormap('jet'),colorbar,set(gca,'FontSize',16);daspect([1 1 1]);
    
    saveas(fig1,char(strcat('B0_invivo_TMSMR28Ch_repro','.tif')));
 %%   
    



