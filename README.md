
# Simulation Studies
Call the simu.R in the simulation folder
```
source('simu.R')
```

# Experiments on Invariance

Run the scripts in the folder invariance


# Data Application
## Download Data
Follow the instructions in http://adni.loni.usc.edu/ to download the image data. A list of images and subjects used in the paper [Intrinsic Riemannian Functional Data Analysis for Sparse Longitudinal Observations](https://arxiv.org/abs/2009.07427) is given in the file subject-list.txt.

## Preparation

Install cmake:  
```
sudo apt install cmake

```

Install dcm2niix


```
cd ~
git clone https://github.com/rordenlab/dcm2niix.git
cd dcm2niix
mkdir build && cd build
cmake -DZLIB_IMPLEMENTATION=Cloudflare -DUSE_JPEGLS=ON -DUSE_OPENJPEG=ON ..
make
```

Install conda: Follow the instruction here: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

Install python2

```
sudo apt install python

```

Install FSL: First, follow the instructions of https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation to download ``fslinstaller.py’’ into ~/tmp. Then, execute the following commands

```
python2 ~/tmp/fslinstaller.py
Install MRtrix3 (mrcalc, dwinormalise)
conda install -c omnia eigen3
conda install -c dsdale24 qt5
conda install -c mrtrix3 mrtrix3
```

Install ANTs: Follow the instruction: http://neuro.debian.net/install_pkg.html?p=ants


Use Python 3.6 (for the first time running)
```
conda create -n py36 python=3.6 anaconda
sudo apt install python3-pip
```

Install Scilpy: Download scilpy: https://github.com/scilus/scilpy and unzip it into ~/scilpy

```
cd ~/scilpy
conda activate py36 (if py36 has not been activated yet)
pip install -r requirements.txt
pip install -e .
conda install python.app
```

Patching: copy private_mask_adjust.py and scil_apply_bias_field_on_dwiv1.py to ~/scilpy/scripts


## Processing

Suppose all DCM images of a single study of a single subject is in the folder “~/tmp/data”

First, create some working directories:
```
cd ~/tmp
mkdir diffusion/
mkdir diffusion/dti
```

Extract necessary information from DCM:
```
~/dcm2niix/build/bin/dcm2niix -o ~/tmp/diffusion ~/tmp/data
```

Move .nii  b0 bval bvec to diffusion/
```
cd ~/tmp/diffusion
mv *.bval bval
mv *.bvec bvec
mv *.nii temp.nii
```

If py36 has not been activated yet
```
conda activate py36
```

Start processing:
```
python ~/scilpy/scripts/scil_extract_b0.py temp.nii bval bvec b0.nii --mean --b0_thr 10

bet b0.nii b0_bet.nii -m -R -f 0.10

mrcalc temp.nii b0_bet_mask.nii.gz  -mult dwi_bet.nii.gz -quiet -force

N4BiasFieldCorrection -i b0_bet.nii.gz -x b0_bet_mask.nii.gz  -o b0_n4.nii.gz, bias_field_b0.nii.gz -c 300x150x75x50, 1e-6 -v 1

python ~/scilpy/scripts/private_mask_adjust.py bias_field_b0.nii.gz b0_bet_mask_adjusted.nii.gz --mask b0_bet_mask.nii.gz -f

python ~/scilpy/scripts/scil_apply_bias_field_on_dwiv1.py dwi_bet.nii.gz bias_field_b0.nii.gz dwi_n4.nii.gz --mask b0_bet_mask_adjusted.nii.gz -f

python ~/scilpy/scripts/scil_crop_volume.py dwi_n4.nii.gz dwi_cropped.nii.gz -f --output_bbox dwi_boundingBox.pkl -f

python ~/scilpy/scripts/scil_crop_volume.py b0_bet.nii.gz b0_cropped.nii.gz  --input_bbox dwi_boundingBox.pkl -f

python ~/scilpy/scripts/scil_crop_volume.py b0_bet_mask_adjusted.nii.gz b0_mask_cropped.nii.gz        --input_bbox dwi_boundingBox.pkl -f

dwinormalise individual dwi_cropped.nii.gz b0_cropped.nii.gz dwi_normalized.nii.gz -fslgrad bvec bval

python ~/scilpy/scripts/scil_resample_volume.py dwi_normalized.nii.gz dwi_resample.nii.gz             --resolution 1 --interp "cubic" -f

fslmaths dwi_resample.nii.gz -thr 0 dwi_resample_clipped.nii.gz

python ~/scilpy/scripts/scil_resample_volume.py b0_mask_cropped.nii.gz \
            mask_resample.nii.gz  --ref dwi_resample.nii.gz --enforce_dimensions --interp nn -f

mrcalc dwi_resample_clipped.nii.gz mask_resample.nii.gz -mult dwi_resampled.nii.gz -quiet -force

python ~/scilpy/scripts/scil_extract_b0.py dwi_resampled.nii.gz bval bvec b0_resampled.nii.gz --mean --b0_thr 10

mrthreshold b0_resampled.nii.gz b0_mask_resampled.nii.gz --abs 0.00001 -force

python ~/scilpy/scripts/scil_extract_dwi_shell.py dwi_resampled.nii.gz bval bvec 0 1000 dwi_dti.nii.gz bval_dti bvec_dti -t 20 -f

python ~/scilpy/scripts/scil_extract_dwi_shell.py dwi_resampled.nii.gz bval bvec 0 1000 2000 3000 dwi_fodf.nii.gz bval_fodf bvec_fodf -t 20 -f

python ~/scilpy/scripts/scil_compute_dti_metrics.py dwi_dti.nii.gz  bval_dti bvec_dti --mask b0_mask_resampled.nii.gz \
        --not_all \
        --ad dti/ad.nii.gz --evecs dti/evecs.nii.gz\
        --evals dti/evals.nii.gz --fa dti/fa.nii.gz\
        --ga dti/ga.nii.gz --rgb dti/rgb.nii.gz\
        --md dti/md.nii.gz --mode dti/mode.nii.gz\
        --norm dti/norm.nii.gz --rd dti/rd.nii.gz\
        --tensor dti/tensor.nii.gz\
        --non-physical dti/nonphysical.nii.gz\
        --pulsation dti/pulsation.nii.gz\
        -f
```

## Alternative processing 

```
cd ~/tmp
mkdir diffusion/
mkdir diffusion/dti

~/dcm2niix/build/bin/dcm2niix -o ~/tmp/diffusion ~/tmp/data

cd ~/tmp/diffusion
mv *.bval bval
mv *.bvec bvec
mv *.nii temp.nii
rm *.json


mrconvert -fslgrad bvec bval temp.nii temp.mif -force

dwidenoise temp.mif temp_denoised.mif -noise temp_noise.mif -force

dwifslpreproc temp_denoised.mif temp_preproc.mif -rpe_none -pe_dir j -nocleanup -force

dwiextract temp_preproc.mif - -bzero | mrmath - mean b0.mif -axis 3

mrconvert temp_preproc.mif data.nii.gz -export_grad_mrtrix temp_grad_table -export_grad_fsl bvecs bvals -force

mrconvert b0.mif b0.nii.gz

bet2 b0.nii.gz b0_bet.nii.gz -f 0.2

dwibiascorrect ants data.nii.gz data_bias_corrected.nii.gz -fslgrad bvecs bvals -mask b0_bet.nii.gz -ants.c [300x150x75x50,1e-6] -force

dwinormalise individual data_bias_corrected.nii.gz b0_bet.nii.gz data_normalized.nii.gz -fslgrad bvecs bvals

dtifit --data=data_normalized.nii.gz --out=dti --mask=b0_bet.nii.gz --bvecs=bvecs --bvals=bvals --save_tensor

```

Registration:

```
ImageMath 4 dtiComp.nii.gz TimeSeriesDisassemble dti_tensor.nii.gz
i=0 
for index in xx xy xz yy yz zz; do
   mv dtiComp100${i}.nii.gz dtiComp_${index}.nii.gz
   i=$((i+1))
done
ImageMath 3 dtAnts.nii.gz ComponentTo3DTensor dtiComp_
ImageMath 3 dti_FA.nii.gz TensorFA dtAnts.nii.gz
```

## Calculate the diffusion tensors
Follow the instruction of hippo.m to derive a difussion tensor for each image. The script produces a MATLAB data file hippo.mat that contains the tensors and the associated image IDs. Now we assume that the meta information of all images is stored in a data frame called 'meta' with four columns: a column of 'image_id', a column of 'subject_id', a column of 'age', and a column of 'group' (CN or AD). Use the following R commands to combine the tensors with metadata, after loading the meta data frame into R:

Then use the following R commands to convert the MATLAB data file into a R data file.
```
library(R.matlab)
hippo.mat <- readMat('hippo.mat')
idx <- match(hippo.mat$imageIds,meta$image_id)
dti=hippo.mat$dti.result
metatinfo=meta[idx,]
save(file='hippo.RData', dti, metainfo)
```

## iRFDA analysis
Now call the analysis.R in the data-application folder:
```
source('analysis-script.R')
```
This will create two R data files, AZ.RData and CN.RData (which can be also found in the data-application folder). Use the script plot.R (depending on 3d.par.2.RData) to generate the plots in the paper.

```
source('plot.R')
```
