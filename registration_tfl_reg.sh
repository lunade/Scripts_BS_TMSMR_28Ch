# Lucia Navarro de Lara, PhD
# MGH, Martinos Center
# 04/25/2023
#lnavarrodelara@mgh.harvard.edu

flirt -in snr_20Ch.nii -ref snr_with_tms.nii -o snr_20Ch_reg.nii -cost mutualinfo -omat 20Ch2TMS.mat
flirt -in tfl_20Ch_reg_gre.nii.gz -ref tfl_with_tms_gre_reg.nii.gz -omat 20Ch2TMS_gre.mat -bins 256 -cost normcorr -searchrx -20 20 -searchry -20 -20 -searchrz -20 20 -dof 6 -interp nearestneighbour -out tfl_20Ch_gre_reg2.nii
flirt -in tfl_20Ch_reg.nii.gz -ref tfl_with_tms_reg.nii.gz -init 20Ch2TMS_gre.mat -applyxfm -interp nearestneighbour -out tfl_20Ch_reg2.nii
