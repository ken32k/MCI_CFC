#!/bin/bash
###
 # @Author: Ken32g ken32k@163.com
 # @Date: 2024-04-09 15:22:23
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-07-01 14:21:38
 # @FilePath: /csvd-sfc/svd_file_ls-subjects.sh
 # @Description: Statistics for middle files in preprocessing
### 

# Function to echo with colored background
function echo_with_color() {
    local color=$1
    shift
    echo -e "\e[48;5;${color}m$@\e[0m"
}

# Path to the folder to scan
mica_out_path=${CSVD_DIR}/csvd_mica_out/micapipe_v0.2.0
bids_path=${CSVD_DIR}/csvd_mica_bids

subs=""
N=0
M=0
# Loop through subfolders
for subfolder in "$bids_path"/*/; do
    sub_id=$(basename "$subfolder")
    # tar_file=(${mica_out_path}/${sub_id}/dwi/*/sub-*_space-dwi_atlas-schaefer-200_desc-iFOD2-2M-SIFT2_full-assignments.txt)
    # tar_file=(${bids_path}/${sub_id}/anat/sub-*_T1w.nii.gz)
    # tar_file=(${bids_path}/${sub_id}/anat/sub-*_FLAIR.nii.gz)
    tar_file=(${mica_out_path}/${sub_id}/func/desc-se_dir-AP_rest_bold/surf/${sub_id}_surf-fsLR-32k_atlas-schaefer-200_desc-FC.shape.gii)
    tar_file=(${bids_path}/${sub_id}/func/sub-*_rest*.nii.gz)
    # tar_file=("$subfolder"/anat/*.nii.gz)
    # Check if subfolder is empty
    if [ ! "$(ls -A "${subfolder}")" ]; then
        # Store subfolder name with yellow background
        subs+="$(echo_with_color 3 "$sub_id") "
    elif [ -f "${tar_file[0]}" ]; then
        # Store subfolder name with blue background
        N=$((N+1))
        subs+="$(echo_with_color 4 "$sub_id") "
    elif [ ! -d "${mica_out_path}/${sub_id}" ]; then
        # Store subfolder name with blue background
        subs+="$(echo_with_color 1 "$sub_id") "
        # M=$((M+1))
    else
        # Store subfolder name with red background
        subs+="$(echo_with_color 2 "$sub_id") "
        M=$((M+1))
    fi
done

# Print all subfolder names in a single line
echo "$subs"

echo ${N}, ${M}