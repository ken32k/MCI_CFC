#!/bin/bash
# Dicom to Nifti conversion for MCI-CFC dataset using svd_file-dcm2nii.sh script. This script will submit a batch job for each subject that does not have the output file already present. 

out=${MCI_CFC_DIR}/mica_out
subjects=(${MCI_CFC_DIR}/mica_bids/sub-*)

batch_subjects=(${subjects[@]:0:100})

N=0
# dest_file=${out}/micapipe_v0.2.0/sub-${sub_id: -3}/dwi/connectomes-wmh-2M/sub-${sub_id: -3}_space-dwi_atlas-schaefer-200_desc-iFOD2-2M-SIFT2_full-assignments.txt
# dest_file="~/absolutelynotexist"
# dest_file="${MCI_CFC_DIR}/mica_bids/sub-${sub_id: -3}/func/sub-${sub_id: -3}_dir-AP_rest_bold.nii.gz"

for subject in "${batch_subjects[@]:0:100}"; do
    sub_id=${subject: -3}
    dest_file="${MCI_CFC_DIR}/mica_bids/sub-${sub_id}/func/sub-${sub_id}_dir-AP_rest_bold.nii.gz"

    if [ "$N" -gt 100 ]; then
        exit
    elif [ ! -f $dest_file ]; then
        N=$((N + 1))
        echo "${sub_id}"

        sbatch --output="${MCI_CFC_DIR}/log/sub-${sub_id}-%j.txt" \
            --error="${MCI_CFC_DIR}/log/sub-${sub_id}-%j.err" \
            --job-name="svd-${sub_id}" \
            --nodes=1 \
            --ntasks=1 \
            --partition=bme_cpu \
            --cpus-per-task=1 \
            --time=04:00:00 \
            ${HOME}/communication//svd_file-dcm2nii.sh "${sub_id}"
    else
        echo "output for ${subject} exist!"
    fi
done
