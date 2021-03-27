#!/bin/bash

export PATH="/home/ulysses/dcm2niix/build/bin:$PATH"

dcmdir=$1
eddy=$2
template=$3
scriptdir=$4

dcm2niix -o . ${dcmdir}

mv *.bval bval
mv *.bvec bvec
mv *.nii temp.nii
rm *.json

mrconvert -fslgrad bvec bval ${eddy} temp.mif -force

dwidenoise temp.mif temp_denoised.mif -noise temp_noise.mif -force

dwiextract temp_denoised.mif - -bzero | mrmath - mean b0.mif -axis 3

mrconvert temp_denoised.mif data.nii.gz -export_grad_mrtrix temp_grad_table -export_grad_fsl bvecs bvals -force

mrconvert b0.mif b0.nii.gz

bet2 b0.nii.gz b0_bet.nii.gz -f 0.2

dwibiascorrect ants data.nii.gz data_bias_corrected.nii.gz -fslgrad bvecs bvals -mask b0_bet.nii.gz -ants.c [300x150x75x50,1e-6] -force

dwinormalise individual data_bias_corrected.nii.gz b0_bet.nii.gz data_normalized.nii.gz -fslgrad bvecs bvals

dtifit --data=data_normalized.nii.gz --out=dti --mask=b0_bet.nii.gz --bvecs=bvecs --bvals=bvals --save_tensor

# registration

ImageMath 4 dtiComp.nii.gz TimeSeriesDisassemble dti_tensor.nii.gz
i=0 
for index in xx xy xz yy yz zz; do
   mv dtiComp100${i}.nii.gz dtiComp_${index}.nii.gz
   i=$((i+1))
done
ImageMath 3 dtAnts.nii.gz ComponentTo3DTensor dtiComp_
ImageMath 3 dti_FA.nii.gz TensorFA dtAnts.nii.gz

${scriptdir}regHelper.sh dti_FA.nii.gz dtAnts.nii.gz ${template} ./

ImageMath 3 aligned_FA.nii.gz TensorFA DTDeformed.nii.gz
