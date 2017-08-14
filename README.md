fMRIANTs
========

If you want a full pipeline, see [ANTsR fMRI example](http://htmlpreview.github.io/?https://github.com/stnava/fMRIANTs/blob/master/ANTsRfMRI_FAQ.html).  Otherwise, see below.

Minimal fMRI pre-processing with ANTS

Let's first define some data and create an average (target) image.
```
fmri=data/bold.nii.gz
nm=data/AFFINE
ref=${nm}_avg.nii.gz
antsMotionCorr -d 3 -a $fmri -o $ref
```

The variable `$ref` is the target image (fixed image) to which we will motion
correct.

Let's do affine motion correction first.
```
antsMotionCorr  -d 3 \
-o [ ${nm}, ${nm}.nii.gz,${nm}_avg.nii.gz] \
-m MI[${ref}, ${fmri}, 1 , 32 , Regular, 0.1  ] \
-t Affine[ 0.1 ] -u 1 -e 1 -s 1x0 -f 2x1 \
-i 15x3 -n 3  -w 1
```
This is a 'fast' example - you might change `-i` parameters
to something larger and `Regular, 0.1` to `Regular, 0.2` for 'real' data.

The parameter `-w 1` means write out the displacement field.  This will
write a 4D displacement field that captures the affine induced motion
at each voxel.  This can be used to concatenate space-time transformations
or to control for motion, statistically, in a spatially varying way.

Now we create a 4D target image for 4D motion deformable motion correction.
First, parse the header information to find out the number of time points.
```
hislice=`PrintHeader $fmri | grep Dimens | cut -d ',' -f 4 | cut -d ']' -f 1`
tr=`PrintHeader $fmri | grep "Voxel Spac" | cut -d ',' -f 4 | cut -d ']' -f 1`
```

Replicate the fixed 3D image `hislice` times to make a new 4D fixed image.
```
fxd=${nm}_fixed.nii.gz
ImageMath 3 $fxd ReplicateImage ${nm}_avg.nii.gz $hislice $tr 0
```

Finally, run SyN to map the time series to this fixed space.
```
tx=SyN[0.1,3,0.0] # critical parameters (though others matter too)
antsRegistration --dimensionality 4 -r ${nm}Warp.nii.gz \
      --output   [${nm}_,${nm}Warped.nii.gz] \
      --interpolation Linear --use-histogram-matching 1 \
      --winsorize-image-intensities [0.005,0.995] --transform $tx \
      --metric meansquares[${fxd},$fmri,1] \
      --convergence [15x2,1e-6,4] --shrink-factors 2x1 \
      --smoothing-sigmas 1x0vox --restrict-deformation 1x1x1x0
```
The parameter `--restrict-deformation 1x1x1x0` prevents deformation across time.

Finally, we generate a 3D transformation to a template then replicate these maps s.t. they can be applied to the original 4D dataset.  This is a useful
example of how ANTs works, in general, to minimize the number of interpolations
that one must apply to the data.
```
antsRegistrationSyNQuick.sh -d 3 -f data/template.nii.gz -m ${nm}_avg.nii.gz \
  -o ${nm}_diff -t s
```

Use `antsApplyTransforms` to combine the displacement field and affine matrix
into a single concatenated transformation stored as a displacement field.
```
# collapse the transformations to a displacement field
antsApplyTransforms -d 3 -o [${nm}_diffCollapsedWarp.nii.gz,1] \
  -t ${nm}_diff1Warp.nii.gz -t ${nm}_diff0GenericAffine.mat \
  -r data/template.nii.gz
```

Replicate the 3D template to 4D and apply all the transformations to the
original BOLD data.
```
ImageMath 3 ${nm}_diff4DCollapsedWarp.nii.gz ReplicateDisplacement \
  ${nm}_diffCollapsedWarp.nii.gz $hislice $tr 0
ImageMath 3 data/template_replicated.nii.gz ReplicateImage \
  data/template.nii.gz $hislice $tr 0
# apply to original bold
antsApplyTransforms -d 4 -o ${nm}_bold2template.nii.gz \
  -t ${nm}_diff4DCollapsedWarp.nii.gz -t ${nm}_0Warp.nii.gz  \
  -r data/template_replicated.nii.gz -i ${nm}Warped.nii.gz
```

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
```

For use in context with T1 processing : see

[NeuroBattery](http://jeffduda.github.io/NeuroBattery/)
