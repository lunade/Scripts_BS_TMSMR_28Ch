

% SNR and g-factor calculation for 2 coils with or without normalization of the FA maps 
% Lucia Navarro de Lara, PhD
% MGH, Martinos Center
% 04/25/2023
%
% Instructions:
% DATA FOLDER SHOULD HAVE 
% For DICOM DATA:
% /DICOM/TFL_XX (FA maps), /DICOM/TFL_GRE_XX (TFL GRE) and /DICOM/GRE_XX ( for snr)
% For Raw Data:
%/rawdata/gre.dat (raw data of a gre image) and /rawdata/noise.dat (noise
% data)

%PROCESSING FOLDER SHOULD have an empty /nifti folder
 
clear all, close all

%Define slices to plot
sli_to_plot = 5:4:60;

filepath = fileparts(mfilename('fullpath'));

% Base folder where data is stored
%base_folder = '/autofs/space/tmsbox_002/users/lucia/SW2release/Data_BS_paper/SNR_invivo'
base_folder = '/DATA_FOLDER';

%Base porcessing folder to save post-processing data
%base_proce_folder = '/autofs/space/tmsbox_002/users/lucia/Processing_TMSMR_28Ch/SNR_comp_April_2023/'
base_proce_folder = '/PROCESSING_FOLDER'

% You must have a nifti order in the processing path, where converted nii
% will be stored. This will be cleared from the beginning
delete([base_proce_folder,'nifti/*']);

%% TMSRM 28Ch with TMS ===============================================================
 % specify path to twix data (raw data from Siemens scanners) for signal and noise 

sig_path_with_tms= [base_folder,filesep,'rawdata',filesep,'gre_TMSMR.dat'];
noise_path_with_tms = [base_folder,filesep,'rawdata',filesep,'gre_TMSMR_noise.dat'];
% 
[snr_cov_with_tms, snr_rss_with_tms, img_with_tms, noisecov_with_tms, noisecorr_with_tms, kspace] = snr_maps_from_raw_data_3D_jrp(sig_path_with_tms, noise_path_with_tms,0);
 
 % it is interesting to compare the optimal covariance-matrix based coil
 % combination with the root-sum-squares combination.  The 'cov'
 % combination uses the noise covariance matrix to weight noisy channels
 % less and remove noise correlations.  High channel-count arrays may
 % benefit from the optimal coil combination
 snr_cov_with_tms = (snr_cov_with_tms);     snr_rss_with_tms = (snr_rss_with_tms);

 % specify path to b1 maps (tfl) and to the DICOMs for the gradient echo scan
 % used to acquire the SNR maps (DICOMs are only used for registration)

tfl_with_tms_filepath = [base_folder,'/DICOM/TFL_B1MAP_45SL_CORONAL_2X2X3_0012'];
tfl_with_tms_gre_filepath = [base_folder,'/DICOM/TFL_B1MAP_45SL_CORONAL_2X2X3_0011'];
snr_with_tms_filepath = [base_folder,'/DICOM/GRE_COIL_MAP_2ISO_FULLHEAD_TR50MS_FA90_0009'];

 

%% NO TMS 32 Ch===================================
% specify path to twix data for signal and noise 

sig_path_32Ch = [base_folder,filesep,'rawdata',filesep,'gre_32Ch.dat'];
noise_path_32Ch  = [base_folder,filesep,'rawdata',filesep,'gre_noise_32Ch.dat'];
[snr_cov_32Ch , snr_rss_32Ch , img_32Ch , noisecov_32Ch , noisecorr_32Ch , kspace ] = snr_maps_from_raw_data_3D_jrp(sig_path_32Ch , noise_path_32Ch ,0);
 snr_cov_32Ch  = (snr_cov_32Ch );     snr_rss_32Ch  = (snr_rss_32Ch );

 % specify path to b1 maps (tfl) and to the DICOMs for the gradient echo scan
 % used to acquire the SNR maps (DICOMs are only used for registration)


tfl_32Ch_filepath = [base_folder,'/DICOM/TFL_B1MAP_45SL_CORONAL_2X2X3_0020'];
tfl_32Ch_gre_filepath = [base_folder,'/DICOM/TFL_B1MAP_45SL_CORONAL_2X2X3_0019'];
snr_32Ch_filepath = [base_folder,'/DICOM/GRE_COIL_MAP_2ISO_FULLHEAD_TR50MS_FA90_0017'];
%% convert DICOMs to NIFTI format =========================================



dcm2niix_wrapper_l2(tfl_with_tms_filepath,[base_proce_folder,'nifti'],'tfl_with_tms')
dcm2niix_wrapper_l2(tfl_with_tms_gre_filepath,[base_proce_folder,'nifti'],'tfl_with_tms_gre')
dcm2niix_wrapper_l2(tfl_32Ch_filepath,[base_proce_folder,'nifti'],'tfl_32Ch')
dcm2niix_wrapper_l2(tfl_32Ch_gre_filepath,[base_proce_folder,'nifti'],'tfl_32Ch_gre')
dcm2niix_wrapper_l2(snr_with_tms_filepath,[base_proce_folder,'nifti'],'gre_with_tms');
dcm2niix_wrapper_l2(snr_32Ch_filepath,[base_proce_folder,'nifti'],'gre_32Ch');        


%% replace original gre with calculated snr niftis

unix(['cp ',[base_proce_folder,'nifti/gre_with_tms.nii'],' ',[base_proce_folder,'nifti/snr_with_tms.nii']]),
% snr has to be flipped to have same orientation as original data
% the snr eval results are bigger than original data. The extra slices
% have to be cropped out to avoid issues later on
snr_cov_with_tms_f=fliplr(snr_cov_with_tms(:,:,5:end-4));
nifti_cheat_replace_vol([base_proce_folder,'nifti/snr_with_tms.nii'],squeeze(rot90(snr_cov_with_tms_f,-1)));
unix(['cp ',[base_proce_folder,'nifti/gre_32Ch.nii'],' ',[base_proce_folder,'nifti/snr_32Ch.nii']]),
% snr has to be flipped to have same orientation as original data
snr_cov_32Ch_f=fliplr(snr_cov_32Ch(:,:,5:end-4));
nifti_cheat_replace_vol([base_proce_folder,'nifti/snr_32Ch.nii'],squeeze(rot90(snr_cov_32Ch_f,-1)));
  
% make copies of niftis for later use in overwriting volumes
unix(['cp ',[base_proce_folder,'nifti/snr_with_tms.nii'],' ',[base_proce_folder,'nifti/snr_with_tms_rss.nii']]), 
unix(['cp ',[base_proce_folder,'nifti/snr_32Ch.nii'],' ', [base_proce_folder,'nifti/snr_32Ch_rss.nii']]);
snr_rss_with_tms_f=fliplr(snr_rss_with_tms(:,:,5:end-4));
snr_rss_32Ch_f=fliplr(snr_rss_32Ch(:,:,5:end-4));
nifti_cheat_replace_vol([base_proce_folder,'nifti/snr_32Ch_rss.nii'],squeeze(rot90(snr_rss_32Ch_f,-1)));

nifti_cheat_replace_vol([base_proce_folder,'nifti/snr_with_tms_rss.nii'],squeeze(rot90(snr_rss_with_tms_f,-1)));

%% interpolate flip angle maps into coordinate system of snr maps for 32Ch

%run outside on a shell... problem with format of produced matrix
%use commands in 
% registration_tfl.sh
% 
%unix(['flirt -in ',base_proce_folder,'nifti/tfl_with_tms_gre.nii -ref ',base_proce_folder,'nifti/gre_with_tms.nii -omat ',base_proce_folder,...
%     'nifti/tfl2gre_2.mat -bins 256 -cost normcorr -searchrx -20 20 -searchry -20 -20 -searchrz -20 20 -dof 6 -interp nearestneighbour -out ',base_proce_folder,'nifti/tfl_with_tms_gre_reg2.nii']),
% unix(['flirt -in ',base_proce_folder,'nifti/tfl_with_tms.nii -ref ',base_proce_folder,'nifti/gre_with_tms.nii -init ',base_proce_folder,...
%     'nifti/tfl2gre_2.mat -applyxfm -interp nearestneighbour -out ',base_proce_folder,'nifti/tfl_with_tms_reg2.nii']),
%%

tfl_with_tms_reg = (load_nifti([base_proce_folder,'nifti/tfl_with_tms_reg.nii.gz']));  tfl_with_tms_reg.vol=(tfl_with_tms_reg.vol);

tfl_32Ch_reg = (load_nifti([base_proce_folder,'nifti/tfl_32Ch_reg.nii.gz']));  tfl_32Ch_reg.vol=(tfl_32Ch_reg.vol);
%% Normalization step  
% No flip angle normalization in this case. 
snr_with_tms_norm = snr_cov_with_tms_f;
snr_32Ch_norm = snr_cov_32Ch_f;
%If NORMALIZATION is required,
% UNCOMMENT THESE LINES!!!  
% snr_with_tms_norm = ((snr_cov_with_tms_f./sin(pi*rot90(tfl_with_tms_reg.vol)/10/180)));
% snr_32Ch_norm = ((snr_cov_32Ch_f./sin(pi*rot90(tfl_32Ch_reg.vol)/10/180)));


snr_with_tms_norm = clean_snr(snr_with_tms_norm);
snr_32Ch_norm = clean_snr(snr_32Ch_norm);
nifti_cheat_replace_vol([base_proce_folder,'nifti/snr_with_tms.nii'],squeeze(rot90(snr_with_tms_norm,-1)));
nifti_cheat_replace_vol([base_proce_folder,'nifti/snr_32Ch.nii'],squeeze(rot90(snr_32Ch_norm,-1)));
  
%%% SNR maps are not registered until here  
  figure(401),imagesc([vol2mos(snr_with_tms_norm) vol2mos(snr_32Ch_norm)],[0 200]),colormap(jet),axis image off
%% register SNR maps for the 32Ch array onto reference case TMSMR 20Ch
%run outside on a shell... problem with format of produced matrix
%use registration_tfl_reg.sh
%unix(['flirt -in ',base_proce_folder,'nifti/snr_32Ch.nii -ref ',base_proce_folder,'nifti/snr_with_tms.nii -o ',base_proce_folder,'nifti/snr_32Ch_reg.nii -cost mutualinfo -v -datatype double'])
%%
temp = (load_nifti([base_proce_folder,'nifti/snr_32Ch_reg.nii.gz']));snr_32Ch_norm_final = rot90(temp.vol,1);%figure(32),imagesc(vol2mos(snr_no_tms_norm_final),[0 1000])
tfl_32Ch_reg_final = (load_nifti([base_proce_folder,'nifti/tfl_32Ch_reg2.nii.gz']));  tfl_32Ch_reg_final=rot90(tfl_32Ch_reg_final.vol,1);
%snr_with_tms_norm_final = (snr_with_tms_norm);
% registere
figure(402),imagesc([vol2mos(snr_with_tms_norm) vol2mos(snr_32Ch_norm_final)],[0 200]),colormap(jet),axis image off
 
%%  
% do g-factor computation
% % crop SNR maps to get a tight field of view
crop1 = 20;  crop2 = 1;
img_gfactor_with_tms = img_with_tms(crop2:end-crop2,crop1:end-crop1,:,:);
img_gfactor_32Ch = img_32Ch(crop2:end-crop2,crop1:end-crop1,:,:);

sli=round(34);  
input_crop = 1;
offset=4;

% im_32 = convert_jpr_format_into_multichannel_image(raw_temp_32);
[gfactor_with_tms gfactor_mean_with_tms]=calculate_g_factor((img_gfactor_with_tms(input_crop:end-input_crop,input_crop:end-input_crop+1,:,:)), noisecov_with_tms, sli); 



% im_no_tms = convert_jpr_format_into_multichannel_image(raw_temp_no_tms);
[gfactor_32Ch gfactor_mean_32Ch]=calculate_g_factor((img_gfactor_32Ch(input_crop:end-input_crop,input_crop:end-input_crop+1,:,:)), noisecov_32Ch, sli+offset);


%
figure(103),imagesc([vol2mos((gfactor_with_tms(:,1:1532))); ],[0,1]),colorbar,colormap(jet),axis image off,h=colorbar,set(h,'FontSize',24),daspect([22 10 1]),title('G-factor TMSMR 28CH');

figure(104),imagesc([vol2mos((gfactor_32Ch(:,1:1520))); ],[0,1]),colorbar,colormap(jet),axis image off,h=colorbar,set(h,'FontSize',24),daspect([22 10 1]),title('G-factor 32Ch');
%%
figure(105),imagesc([vol2mos((gfactor_with_tms(:,1520:2520))); ],[0,1]),colorbar,colormap(jet),axis image off,h=colorbar,set(h,'FontSize',24),daspect([22 10 1]),title('G-factor TMSMR 28CH');

figure(106),imagesc([vol2mos((gfactor_32Ch(:,1520:2520))); ],[0,1]),colorbar,colormap(jet),axis image off,h=colorbar,set(h,'FontSize',24),daspect([22 10 1]),title('G-factor 32Ch');



 
%%  SAVE RESULTS (optional)
 
save([base_proce_folder,'/snr_mapping_TMSMR_32Ch_comp_April_2023.mat'],'snr_cov_32Ch','snr_32Ch_norm','snr_32Ch_norm_final',...
    'snr_with_tms_norm','snr_rss_with_tms','snr_cov_with_tms','img_with_tms','img_32Ch','noisecov_32Ch', 'noisecorr_32Ch', 'noisecov_with_tms', 'noisecorr_with_tms')
save([base_proce_folder,'/field_mapping_TMSMR_32Ch_comp_April_2023.mat'], 'tfl_with_tms_reg', 'tfl_32Ch_reg','tfl_32Ch_reg_final');
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'gfactor_with_tms', 'gfactor_32Ch', 'gfactor_mean_with_tms',  'gfactor_mean_32Ch')  
return;

   %% =====================================================================
   % supporting functions
   
function [scaled_snr_map, snr_map_interp, tfl, sin_tfl]  = handle_interp_scale_snr(path1, path2)

eval(['unix(''flirt -in ./nifti/',path1,'.nii -ref ./nifti/',path2,'.nii -omat ./nifti/invol2refvol.mat -interp nearestneighbour -applyxfm -usesqform -out ./nifti/interp_out.nii -v'')']);
   snr_map_interp=load_nifti('./nifti/interp_out.nii.gz');  snr_map_interp=snr_map_interp.vol;
   
   tfl = load_nifti(['./nifti/',path2,'.nii']);
   tfl=tfl.vol;
   sin_tfl = sin(tfl/10*pi/180);
   sin_tfl(sin_tfl==0)=1;
   scaled_snr_map = snr_map_interp./sin_tfl;

end




function out=interp_pe(in)
    temp=in;
        for aa=1:numel(temp(1,:,1))
            for bb=1:numel(temp(1,1,:))
                out(:,aa,bb) = interp(in(:,aa,bb),2);
            end
        end

end


   
function [gfactor gfactor_mean]  = calculate_g_factor(im, noise,sli)

im=interp_pe(squeeze(im(:,:,sli,:)));
im=rot90(im);

mag_temp=squeeze(im);  mag_temp = conj(mag_temp).*mag_temp;
mag=sum(mag_temp,3);
mask=mag>.3*mean(mean(mag));  mask=double(mask);
mask_zeros=logical(mask);
mask(mask==0)=nan;

% gfactor=0;

gfactor1 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,1,2);  gfactor1=fix_size(im(:,:,1),gfactor1); gfactor1=1./gfactor1.*mask;
gfactor2 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,1,3);  gfactor2=fix_size(im(:,:,1),gfactor2);  gfactor2=1./gfactor2.*mask;
gfactor3 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,1,4);  gfactor3=fix_size(im(:,:,1),gfactor3);  gfactor3=1./gfactor3.*mask;
gfactor4 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,1,5);  gfactor4=fix_size(im(:,:,1),gfactor4);  gfactor4=1./gfactor4.*mask;
gfactor5 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,1,6);  gfactor5=fix_size(im(:,:,1),gfactor5);  gfactor5=1./gfactor5.*mask;
gfactor6 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,1,7);  gfactor6=fix_size(im(:,:,1),gfactor6);  gfactor6=1./gfactor6.*mask;

gfactor7 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,2,2);  gfactor7=fix_size(im(:,:,1),gfactor7); gfactor7=1./gfactor7.*mask;
gfactor8 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,3,3);  gfactor8=fix_size(im(:,:,1),gfactor8); gfactor8=1./gfactor8.*mask;
gfactor9 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,4,4);  gfactor9=fix_size(im(:,:,1),gfactor9); gfactor9=1./gfactor9.*mask;
gfactor10 = mrir_array_SENSE_gfactor_2d(squeeze(im),noise,5,5);  gfactor10=fix_size(im(:,:,1),gfactor10); gfactor10=1./gfactor10.*mask;


gfactor = [gfactor1 gfactor2 gfactor3 gfactor4 gfactor5 gfactor6 gfactor7 gfactor8 gfactor9 gfactor10];

temp=gfactor1(mask_zeros(:));  temp(isinf(temp))=0;  gfactor_mean(1) = mean(temp);
temp=gfactor2(mask_zeros(:)); temp(isinf(temp))=0;  gfactor_mean(2) = mean(temp);
temp=gfactor3(mask_zeros(:));  temp(isinf(temp))=0;gfactor_mean(3) = mean(temp);
temp=gfactor4(mask_zeros(:)); temp(isinf(temp))=0; gfactor_mean(4) = mean(temp);
temp=gfactor5(mask_zeros(:));  temp(isinf(temp))=0; gfactor_mean(5) = mean(temp);
temp=gfactor6(mask_zeros(:)); temp(isinf(temp))=0; gfactor_mean(6) = mean(temp);
temp=gfactor7(mask_zeros(:)); temp(isinf(temp))=0; gfactor_mean(7) = mean(temp);
temp=gfactor8(mask_zeros(:)); temp(isinf(temp))=0; gfactor_mean(8) = mean(temp);
temp=gfactor9(mask_zeros(:)); temp(isinf(temp))=0; gfactor_mean(9) = mean(temp);
temp=gfactor10(mask_zeros(:)); temp(isinf(temp))=0; gfactor_mean(10) = mean(temp);


end


function gfactor_out=fix_size(im,gfactor)
dims1=size(im);  dims2=size(gfactor);
dims_diff=dims1-dims2;
gfactor_out=gfactor;

if dims_diff(1)>0
    gfactor_out(end:end + dims_diff(1),:) = 0;
end

if dims_diff(2)>0
    gfactor_out(:, end:end + dims_diff(2)) = 0;
end


end

function out = get_average_noise_corr(in)

    temp = triu(in,1);
    num = sum(sum(abs(temp)>0));
    out=100*sum(abs(temp(:))./num);

end

function add_text

text(150,240,'nova_stx 20Ch','color','white','FontSize',24),
text(650,240,'no_tmsch','color','white','FontSize',24),


end

function out = clean_snr(in)
 
    in(isinf(in)) = 0;
    in(in>5000)=0;
    out=in;

end

