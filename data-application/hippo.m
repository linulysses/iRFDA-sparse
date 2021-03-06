% Probability masks for hippocumpus region
LEFT_HIPPO_PATH = 'template/harvardoxford-subcortical_prob_Left Hippocampus.nii.gz';
RIGHT_HIPPO_PATH = 'template/harvardoxford-subcortical_prob_Right Hippocampus.nii.gz';

% The root of preprocessed data
% Each subfolder (named by the image ID) of this folder contains
% preprocessed images and files generated by preprocessing pipeline
DATA_DIR = ''; 

% Template used, available from FSL
TEMPLATE_PATH = 'FMRIB58_FA_1mm.nii.gz';

% Number of cores for parallel computing
N_CORE = 16;

%% 
lh = niftiread(LEFT_HIPPO_PATH);
rh = niftiread(RIGHT_HIPPO_PATH);

% thresholding

hippo_mask = (lh>=80) | (rh>=80);

% FA template
template = double(niftiread(TEMPLATE_PATH));
template = template / max(max(max(template)));

% get the IDs of the images that have been processed
tmp = dir(DATA_DIR);
imageIds = [];
for i = 3:length(tmp)
    fname = tmp(i).name;
    flist = dir([DATA_DIR fname]);
    if length(flist) == 12
        imageIds = [imageIds str2num(fname)]; %#ok<AGROW,ST2NM>
    end
end

n = length(imageIds);

dti_result = zeros(3,3,n);
fa_result = zeros(1,n);
good_image = false(1,n);

mypool = parpool(N_CORE);

parfor i = 1:n
%for i = 1:n
    
    fa_path = [DATA_DIR num2str(imageIds(i)) filesep 'FA.nii.gz'];
    if isfile(fa_path)
        FA = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'FA.nii.gz']);
        fa_result(i) = sum(sum(sum(abs(FA-template))));
    else
        fa_result(i) = Inf;
    end
    
    if fa_result(i) < 2.5 * 1e5 % images with overall FA exceeding this threshold are often not properly acquired.
        good_image(i) = true;
    end
    
    if good_image(i)
        V1 = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'V1Deformed.nii.gz']);
        V2 = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'V2Deformed.nii.gz']);
        V3 = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'V3Deformed.nii.gz']);

        L1 = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'L1Deformed.nii.gz']);
        L2 = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'L2Deformed.nii.gz']);
        L3 = niftiread([DATA_DIR num2str(imageIds(i)) filesep 'L3Deformed.nii.gz']);

        M = lc_mean(V1,V2,V3,L1,L2,L3,hippo_mask);
        dti_result(:,:,i) = M;
    else
        dti_result(:,:,i) = Inf;
    end
    
    if rem(i,50) == 0
        disp(['i=' num2str(i) ': completed']);
    end
    
end

dti_result = dti_result(:,:,good_image);
imageIds = imageIds(good_image);
fa_result = fa_result(good_image);

save('hippo.mat','hippo_mask','dti_result','imageIds','fa_result');

delete(mypool);