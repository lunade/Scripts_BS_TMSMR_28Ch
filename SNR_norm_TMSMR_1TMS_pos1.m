% Script to calculate SNR dif
% 
% Lucia Navarro de Lara, PhD
% MGH, Martinos Center
% 04/25/2023
%
% Instructions:
% In this case the SNR commands used are for 2D sequences
% PROCESSING FOLDER folder to save processed data
% DATA_FOLDER is where data is stored
% DATA FOLDER SHOULD HAVE 
% For DICOM DATA:
% /DICOM/TFL_XX (FA maps), /DICOM/TFL_GRE_XX (TFL GRE) and /DICOM/GRE_XX ( for snr)
% For Raw Data:
%/rawdata/gre.dat (raw data of a gre image) and /rawdata/noise.dat (noise
% data) 
% mask_manual2.mat is used to mask the maps
% FIGURE_FOLDER is the folder to save figures

%% TMSMR snr normalized to B1 plus map angle ( to sin (B1))
% Works only with 2D GRE!!!
close all; clear all;
data_folder='DATA_FOLDER'
process_folder= 'PROCESSING_FOLDER'
%%
pathrootdata= [data_folder '/rawdata'];
pathrootdata_B1=[data_folder '/DICOM'];
pathrootprocess = [process_folder];

noise_file_TMSMR_withTMS = [pathrootdata '/gre_withTMS_pos1_noise.dat'];
noise_file_TMSMR_noTMS = [pathrootdata '/gre_noTMS_noise.dat'];
image_file_TMSMR_withTMS = [pathrootdata '/gre_withTMS_pos1.dat'];%%
image_file_TMSMR_withTMS_filepath = [pathrootdata_B1 '/GRE_COIL_MAP_2X2X3_FULLHEAD_45_TMS_1_POS1_0011'];
image_file_TMSMR_noTMS = [pathrootdata '/gre_noTMS.dat'];%%
image_file_TMSMR_noTMS_filepath = [pathrootdata_B1 '/GRE_COIL_MAP_2X2X3_FULLHEAD_45_0007'];
%% Eval withTMS JP call functions
meas1 = read_meas_dat__fast(image_file_TMSMR_withTMS);
meas2 = read_meas_dat__fast(noise_file_TMSMR_withTMS)   % 90 deg flip



[snr_rss_TMSMR_withTMS, snr_cov_TMSMR_withTMS, noise_cov_TMSMR_withTMS, noise_corr_TMSMR_withTMS] = generate_cov_and_snr(meas1, meas2); 
snr_rss_TMSMR_withTMS=squeeze(snr_rss_TMSMR_withTMS); snr_cov_TMSMR_withTMS=squeeze(snr_cov_TMSMR_withTMS);
%% Eval TMSMR no TMS
meas3 = read_meas_dat__fast(image_file_TMSMR_noTMS)   % 90 deg flip
meas4 = read_meas_dat__fast(noise_file_TMSMR_noTMS)   % noise



[snr_rss_TMSMR_noTMS, snr_cov_TMSMR_noTMS, noise_cov_TMSMR_noTMS, noise_corr_TMSMR_noTMS] = generate_cov_and_snr(meas3, meas4); 
snr_rss_TMSMR_noTMS=squeeze(snr_rss_TMSMR_noTMS); snr_cov_TMSMR_noTMS=squeeze(snr_cov_TMSMR_noTMS);


%% Order SNR map Might not be necessary, depending on sequence
snr_cov_noTMS_ordered(:,:,1:2:45)=snr_cov_TMSMR_noTMS(:,:,1:23);snr_cov_noTMS_ordered(:,:,2:2:45)=snr_cov_TMSMR_noTMS(:,:,24:45);
snr_cov_withTMS_ordered(:,:,1:2:45)=snr_cov_TMSMR_withTMS(:,:,1:23);snr_cov_withTMS_ordered(:,:,2:2:45)=snr_cov_TMSMR_withTMS(:,:,24:45);

snr_rss_noTMS_ordered(:,:,1:2:45)=snr_rss_TMSMR_noTMS(:,:,1:23);snr_rss_noTMS_ordered(:,:,2:2:45)=snr_rss_TMSMR_noTMS(:,:,24:45);
snr_rss_withTMS_ordered(:,:,1:2:45)=snr_rss_TMSMR_withTMS(:,:,1:23);snr_rss_withTMS_ordered(:,:,2:2:45)=snr_rss_TMSMR_withTMS(:,:,24:45);

snr_cov_noTMS_ordered_perm=permute(snr_cov_noTMS_ordered,[2,1,3]);
snr_cov_withTMS_ordered_perm=permute(snr_cov_withTMS_ordered,[2,1,3]);
snr_rss_noTMS_ordered_perm=permute(snr_rss_noTMS_ordered, [2,1,3]);
snr_rss_withTMS_ordered_perm=permute(snr_rss_withTMS_ordered, [2,1,3]);
%%
% %% Load B1 plus map with and without TMS
B1plus_withTMS_filepath = [pathrootdata_B1 '/TFL_B1MAP_45SL_CORONAL_2X2X3_TMS1_POS1_0014']; 
B1plus_noTMS_filepath = [pathrootdata_B1 '/TFL_B1MAP_45SL_CORONAL_2X2X3_0010'];
B1plus_withTMS_image_filepath = [pathrootdata_B1 '/TFL_B1MAP_45SL_CORONAL_2X2X3_TMS1_POS1_0013']; 
B1plus_noTMS_image_filepath = [pathrootdata_B1 '/TFL_B1MAP_45SL_CORONAL_2X2X3_0009']; 
%%
dcm2niix_wrapper(B1plus_withTMS_filepath,'nifti','tfl_with_tms_pos1');
dcm2niix_wrapper(B1plus_withTMS_image_filepath,'nifti','tfl_with_tms_pos1_gre');
dcm2niix_wrapper(B1plus_noTMS_filepath,'nifti','tfl_no_tms');
dcm2niix_wrapper(B1plus_noTMS_image_filepath,'nifti','tfl_no_tms_gre');
%%
dcm2niix_wrapper(image_file_TMSMR_withTMS_filepath,'nifti','snr_with_tms_pos1');
dcm2niix_wrapper(image_file_TMSMR_noTMS_filepath,'nifti','snr_no_tms');

%% replace original gre with calculated snr niftis
nifti_cheat_replace_vol([image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1.nii'],squeeze(rot90(snr_cov_withTMS_ordered_perm,2)));
nifti_cheat_replace_vol([image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms.nii'],squeeze(rot90(snr_cov_noTMS_ordered_perm,2)));

% make copies of niftis for later use in overwriting volumes
  unix(['cp ',[image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1.nii'],' ',[image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1_cov.nii']]), 
  unix(['cp ',[image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms.nii'],' ', [image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms_cov.nii']]);

nifti_cheat_replace_vol([image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms.nii'],squeeze(rot90(snr_rss_noTMS_ordered_perm,2)));
nifti_cheat_replace_vol([image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1.nii'],squeeze(rot90(snr_rss_withTMS_ordered_perm,2)));
% make copies of niftis for later use in overwriting volumes
unix(['cp ',[image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1.nii'],' ',[image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1_rss.nii']]), 
unix(['cp ',[image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms.nii'],' ', [image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms_rss.nii']]),
%  
 %% interpolate flip angle maps into coordinate system of snr maps 
 unix(['flirt -in ',B1plus_withTMS_filepath,'/nifti/tfl_with_tms_pos1.nii -ref ',image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1.nii -omat ',image_file_TMSMR_withTMS_filepath,'/nifti/invol2refvol.mat -interp nearestneighbour -applyxfm -usesqform -out ',B1plus_withTMS_filepath,'/nifti/tfl_with_tms_pos1_reg.nii']),
 tfl_with_tms_pos1_reg = (load_nifti([B1plus_withTMS_filepath,'/nifti/tfl_with_tms_pos1_reg.nii.gz']));  tfl_with_tms_pos1_reg.vol=(tfl_with_tms_pos1_reg.vol);
 
 %% interpolate flip angle maps into coordinate system of snr maps for no_tmsch
 unix(['flirt -in ',B1plus_noTMS_filepath,'/nifti/tfl_no_tms.nii -ref ',image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms.nii -omat ',image_file_TMSMR_noTMS_filepath,'/nifti/invol2refvol.mat -interp nearestneighbour -applyxfm -usesqform -out ',B1plus_noTMS_filepath,'/nifti/tfl_no_tms_reg.nii' ]),
 tfl_no_tms_reg = (load_nifti([B1plus_noTMS_filepath,'/nifti/tfl_no_tms_reg.nii.gz']));   tfl_no_tms_reg.vol=(tfl_no_tms_reg.vol);

 %%
%%
snr_cov_flip_withTMS=flip(snr_cov_withTMS_ordered_perm,1);
snr_cov_flip_noTMS=flip(snr_cov_noTMS_ordered_perm,1);
snr_rss_flip_withTMS=flip(snr_rss_withTMS_ordered_perm,1);
snr_rss_flip_noTMS=flip(snr_rss_noTMS_ordered_perm,1);
% Keep L R correct
% Work on shifts COV
snr_cov_flip_withTMS_shift3(:,:,1:44)=snr_cov_flip_withTMS(:,:,2:45);
snr_cov_flip_withTMS_shift3(:,:,45)=snr_cov_flip_withTMS(:,:,1);
snr_cov_flip_noTMS_shift3(:,:,1:44)=snr_cov_flip_noTMS(:,:,2:45);
snr_cov_flip_noTMS_shift3(:,:,45)=snr_cov_flip_noTMS(:,:,1);

snr_cov_withTMS_ordered_perm_shift(2:110,:,:)=snr_cov_flip_withTMS_shift3(1:109,:,:);
snr_cov_withTMS_ordered_perm_shift(1,:,:)=snr_cov_flip_withTMS_shift3(110,:,:);
snr_cov_noTMS_ordered_perm_shift(2:110,:,:)=snr_cov_flip_noTMS_shift3(1:109,:,:);
snr_cov_noTMS_ordered_perm_shift(1,:,:)=snr_cov_flip_noTMS_shift3(110,:,:);

% Work on shifts RSS
snr_rss_flip_withTMS_shift3(:,:,1:44)=snr_rss_flip_withTMS(:,:,2:45);
snr_rss_flip_withTMS_shift3(:,:,45)=snr_rss_flip_withTMS(:,:,1);
snr_rss_flip_noTMS_shift3(:,:,1:44)=snr_rss_flip_noTMS(:,:,2:45);
snr_rss_flip_noTMS_shift3(:,:,45)=snr_rss_flip_noTMS(:,:,1);

snr_rss_withTMS_ordered_perm_shift(2:110,:,:)=snr_rss_flip_withTMS_shift3(1:109,:,:);
snr_rss_withTMS_ordered_perm_shift(1,:,:)=snr_rss_flip_withTMS_shift3(110,:,:);
snr_rss_noTMS_ordered_perm_shift(2:110,:,:)=snr_rss_flip_noTMS_shift3(1:109,:,:);
snr_rss_noTMS_ordered_perm_shift(1,:,:)=snr_rss_flip_noTMS_shift3(110,:,:);
%%
b1_plus_with_tms= tfl_with_tms_pos1_reg.vol;
b1_plus_no_tms= tfl_no_tms_reg.vol;
b1_plus_dif=(b1_plus_with_tms-b1_plus_no_tms)./b1_plus_no_tms;

%%remove flip angle bias from snr maps by dividing by the sine of the flip angle

snr_cov_with_tms_norm = ((snr_cov_withTMS_ordered_perm_shift./sin(pi*b1_plus_with_tms/10/180)));
%   
 snr_cov_no_tms_norm = ((snr_cov_noTMS_ordered_perm_shift./sin(pi*b1_plus_no_tms/10/180)));
 
snr_rss_with_tms_norm = ((snr_rss_withTMS_ordered_perm_shift./sin(pi*b1_plus_with_tms/10/180)));
%   
 snr_rss_no_tms_norm = ((snr_rss_noTMS_ordered_perm_shift./sin(pi*b1_plus_no_tms/10/180))); 
%%

%%
% 
snr_cov_with_tms_norm = clean_snr(snr_cov_with_tms_norm);
% 
snr_cov_no_tms_norm = clean_snr(snr_cov_no_tms_norm);
% 
snr_rss_with_tms_norm = clean_snr(snr_rss_with_tms_norm);
% 
snr_rss_no_tms_norm = clean_snr(snr_rss_no_tms_norm);

%% 
% 

%figure(501),imagesc([vol2mos(snr_cov_with_tms_norm) vol2mos(snr_cov_no_tms_norm)],[100 1500]),colormap(jet),axis image off 
ker=1/8*ones(2,2,2);
B1plus_dif_smooth=b1_plus_dif*100;
B1plus_dif_smooth=convn(B1plus_dif_smooth,ker,'same');
B1plus_dif_smooth(~isfinite(B1plus_dif_smooth))=0;

nifti_cheat_replace_vol([image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1_cov.nii'],squeeze(snr_cov_with_tms_norm));
nifti_cheat_replace_vol([image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms_cov.nii'],squeeze(snr_cov_no_tms_norm));
 
nifti_cheat_replace_vol([image_file_TMSMR_withTMS_filepath,'/nifti/snr_with_tms_pos1_rss.nii'],squeeze(snr_rss_with_tms_norm));
nifti_cheat_replace_vol([image_file_TMSMR_noTMS_filepath,'/nifti/snr_no_tms_rss.nii'],squeeze(snr_rss_no_tms_norm));  
%%% SNR maps normalized not registered until here  
  figure(401),imagesc([vol2mos(snr_rss_with_tms_norm) vol2mos(snr_rss_no_tms_norm)],[100 1500]),colormap(jet),axis image off 
  
  figure(501),imagesc([vol2mos(snr_cov_with_tms_norm) vol2mos(snr_cov_no_tms_norm)],[100 1500]),colormap(jet),axis image off 

  figure(601), imagesc([vol2mos(b1_plus_with_tms/10) vol2mos(b1_plus_no_tms/10)],[0 90]),colormap(jet),axis image off
  %%

%%
 msk=snr_cov_no_tms_norm>100;
 
 
 %% SNR dif
snr_rss_dif_norm= ((snr_rss_with_tms_norm-snr_rss_no_tms_norm)./snr_rss_no_tms_norm)*100;
snr_cov_dif_norm= ((snr_cov_with_tms_norm-snr_cov_no_tms_norm)./snr_cov_no_tms_norm)*100;
snr_rss_dif_norm(~isfinite(snr_rss_dif_norm))=0;% get rid of Inf
snr_cov_dif_norm(~isfinite(snr_cov_dif_norm))=0;% get rid of Inf

% % B1 diff
% B1plus_dif=((dTMS_ordered-dnoTMS_ordered)./dnoTMS_ordered)*100;
% B1plus_dif(~isfinite(B1plus_dif))=0;% get rid of Inf
%% 
%%Save data
save([pathrootprocess,'/results_norm_pos1.mat'],'snr_rss_with_tms_norm','snr_rss_no_tms_norm','snr_cov_with_tms_norm','snr_cov_no_tms_norm','snr_cov_dif_norm','snr_rss_dif_norm','b1_plus_with_tms','b1_plus_no_tms','b1_plus_dif','msk');

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



function out = get_average_noise_corr(in)

    temp = triu(in,1);
    num = sum(sum(abs(temp)>0));
    out=100*sum(abs(temp(:))./num);

end

function add_text

text(150,240,'nova_stx 32ch','color','white','FontSize',24),
text(650,240,'no_tmsch','color','white','FontSize',24),


end

function out = clean_snr(in)
 
    in(isinf(in)) = 0;
    in(in>5000)=0;
    out=in;

end

