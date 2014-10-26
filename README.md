fMRIANTs
========

Minimal fMRI pre-processing with ANTS

```
fmri=data/bold.nii.gz
nm=data/AFFINE
ref=${nm}_avg.nii.gz
antsMotionCorr -d 3 -a $fmri -o $ref
antsMotionCorr  -d 3 -o [ ${nm}, ${nm}.nii.gz,${nm}_avg.nii.gz] -m MI[${ref}, ${fmri}, 1 , 32 , Regular, 0.1  ] -t Affine[ 0.1 ] -u 1 -e 1 -s 1x0 -f 2x1 -i 15x3 -n 3  -w 1
# -w 1 means write out the displacement field
#
# this is a 'fast' example - you should change -i parameters to something larger
# and Regular, 0.01 to Regular, 0.2 for 'real' data -
# see AFFINE_AND_DEFORMABLE below
#
###
### use antsRegistration - must create 4D target
###
hislice=`PrintHeader $fmri | grep Dimens | cut -d ',' -f 4 | cut -d ']' -f 1`
tr=`PrintHeader $fmri | grep "Voxel Spac" | cut -d ',' -f 4 | cut -d ']' -f 1`
fxd=${nm}_fixed.nii.gz
stackavg=" "
for x in `seq 1 $hislice` ; do
  stackavg=" $stackavg ${nm}_avg.nii.gz "
done
ImageMath 4 $fxd TimeSeriesAssemble $tr 0 $stackavg
tx=SyN[0.1,3,0.0] # critical parameters (though others matter too)
antsRegistration --dimensionality 4 -f 1 -r ${nm}Warp.nii.gz \
      --output   [${nm}_,${nm}Warped.nii.gz] \
      --interpolation Linear --use-histogram-matching 1 \
      --winsorize-image-intensities [0.005,0.995] --transform $tx \
      --metric meansquares[${fxd},$fmri,1] \
      --convergence [15x2,1e-6,4] --shrink-factors 2x1 \
      --smoothing-sigmas 1x0vox --restrict-deformation 1x1x1x0
### TODO - show a 3D transformation to a template then pad these maps s.t.
### they can be applied to a 4D dataset - 4D affine and 4D displacement
###########################################
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

[NeuroBattery](http://jeffduda.github.io/NeuroBattery/)
