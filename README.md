fMRIANTs
========

Minimal fMRI pre-processing with ANTS

```
fmri=data/bold.nii.gz
nm=data/AFFINE
ref=${nm}_avg.nii.gz
antsMotionCorr -d 3 -a $fmri -o $ref
antsMotionCorr  -d 3 -o [ ${nm}, ${nm}.nii.gz,${nm}_avg.nii.gz] -m MI[${ref}, ${fmri}, 1 , 32 , Regular, 0.01  ] -t Affine[ 0.1 ] -u 1 -e 1 -s 1x0 -f 2x1 -i 3x0 -n 3  
#
# this is a 'fast' example - you should change -i 3x0 to something larger 
# and Regular, 0.01 to Regular, 0.2 for 'real' data - see AFFINE_AND_DEFORMABLE below
#
###########################################
### now do deformable motion correction ###
###########################################
nm=data/AFFINE_AND_DEFORMABLE
antsMotionCorr  -d 3 -o [ ${nm}, ${nm}.nii.gz,${nm}_avg.nii.gz] -m MI[${ref}, ${fmri}, 1 , 32 , Regular, 0.2  ] -t Affine[ 0.1 ] -u 1 -e 1 -s 1x0 -f 2x1 -i 33x20 -n 3  -m MI[${ref}, ${fmri}, 1 , 32 ] -t GaussianDisplacementField[0.15,3,0.5] -i 33x20 -u 1 -e 1 -s 1x0 -f 4x1 -n 3
```

The csv and result images that come out contain the motion corrected image and the motion nuisance variables - but you need to run the affine version to get the motion variables ( that's a small bug in that the 2nd level registration doesnt maintain the 1st level motion results ) .... It might be interesting to also write out the deformation field in order to quantify distortion but that is a separate issue.

Next you can run the CompCor command to estimate physiological nuisance variables.

```
# 
# CompCor needs a mask --- you can get the mask by (might need to alter the value 200 up or down) 
#
ThresholdImage 3  data/AFFINE_avg.nii.gz data/mask.nii.gz 200 1.e9 
ImageMath 3 data/mask.nii.gz ME data/mask.nii.gz 1 
ImageMath 3 data/mask.nii.gz GetLargestComponent data/mask.nii.gz
ImageMath 3 data/mask.nii.gz MD data/mask.nii.gz 2 
ImageMath 3 data/mask.nii.gz ME data/mask.nii.gz 1 
# you may alternatively run CompCor on the motion-corrected data
ImageMath 4 data/OUT.nii.gz CompCorrAuto data/bold.nii.gz data/mask.nii.gz 6 
#
# the important output will be called  OUT_compcorr.csv 
#
```

For T1 processing : see

(NeuroBattery)[http://jeffduda.github.io/NeuroBattery/]
