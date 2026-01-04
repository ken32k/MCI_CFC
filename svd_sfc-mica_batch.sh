#!/bin/bash
###
 # @Author: Ken32g ken32k@163.com
 # @Date: 2023-12-16 20:52:13
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-07-01 14:19:55
 # @FilePath: /csvd-sfc/svd_sfc-mica_batch.sh
 # @Description: Submit tasks to Slurm
### 

export PATH=$PATH:/home/lingbin/miniconda3/pkgs/mrtrix3-3.0.4-h2bc3f7f_0/bin/
mica_dir=/home/lingbin/ke/communication/mica_script

export FREESURFER_HOME=/usr/local/freesurfer/7.4.1
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/home/lingbin/ke/svd/raw/bidstest_result/fastsurfer

# Define the input folder
input_folder=/home/lingbin/ke/svd/raw/bidstest
output_dir=/home/lingbin/ke/svd/raw/bidstest_result

# Scan the input folder for subfolders and create an array
subfolders=()
while IFS= read -r -d '' subfolder; do
    subfolder_name=$(basename "$subfolder")
    subfolder_name_no_prefix="${subfolder_name#sub-}"
    subfolders+=("$subfolder_name_no_prefix")
done < <(find "$input_folder" -maxdepth 1 -type d -name "sub-*" -print0)

echo mica-batch
# Loop through the subfolders and run the command for each
for subfolder in "${subfolders[@]}"
do  
    if [ ! -f "${output_dir}/micapipe_v0.2.0/${subfolder}/dwi/connectomes/*.gii" ]; then
        # Attempt to run the command for the current subfolder
        echo "Processing subfolder $subfolder..."
        if bash ${mica_dir}/svd_sfc-mica_run.sh "$subfolder"; then
            echo "Command succeeded for subfolder $subfolder"
        else
            echo "Command failed for subfolder $subfolder"
            # Handle the error here, such as retrying the command or exiting the script
        fi
    else
        echo -e '\e[1;34m' "Pass $subfolder"
    fi 
done
