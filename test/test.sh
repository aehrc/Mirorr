#!/bin/bash
mkdir -p output
time ../build/bin/mirorr -m Sample_Male_Pelvis_CT.nii.gz -f Sample_Male_Pelvis_MR.nii.gz -b 3 -c 3 -t rigid  -l output/CT_to_MR-rigid.tfm  --save-moving output/CT_in_MR-rigid.nii.gz  --fresh
time ../build/bin/mirorr -m Sample_Male_Pelvis_CT.nii.gz -f Sample_Male_Pelvis_MR.nii.gz -b 3 -c 3 -t affine -l output/CT_to_MR-affine.tfm --save-moving output/CT_in_MR-affine.nii.gz --fresh

