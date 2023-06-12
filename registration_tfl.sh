# Registration of FA maps to GRE images
#Lucia Navarro de Lara, PhD
#MGH, Martinos Center
# 04/25/2023
#lnavarrodelara@mgh.harvard.edu

flirt -in tfl_with_tms_gre.nii -ref gre_with_tms.nii -omat tfl2gre.mat -bins 256 -cost normcorr -searchrx -20 20 -searchry -20 -20 -searchrz -20 20 -dof 6 -interp nearestneighbour -out tfl_with_tms_gre_reg.nii
flirt -in tfl_with_tms.nii -ref gre_with_tms.nii -init tfl2gre.mat -applyxfm -interp nearestneighbour -out tfl_with_tms_reg.nii
flirt -in tfl_20Ch_gre.nii -ref gre_20Ch.nii -omat tfl2gre_20Ch.mat -bins 256 -cost normcorr -searchrx -20 20 -searchry -20 -20 -searchrz -20 20 -dof 6 -interp nearestneighbour -out tfl_20Ch_reg_gre.nii
flirt -in tfl_20Ch.nii -ref gre_20Ch.nii -init tfl2gre_20Ch.mat -applyxfm -interp nearestneighbour -out tfl_20Ch_reg.nii
