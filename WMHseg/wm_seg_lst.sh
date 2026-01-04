#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-09-17 14:55:52
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-09-26 22:18:46
# @FilePath: /csvd-sfc/WMHseg/mci_wmh_seg_lst.sh
# @Description: WMH segmentation using lst ai
###
sub_num=$1 # 001
export PROJ_HOME=$CSVD_DIR
# Activate LST environment
source ~/apps/lst/bin/activate
# Add Greedy
export PATH=$PATH:~/apps

bids_dir=${PROJ_HOME}/csvd_mica_bids
out_dir=${PROJ_HOME}/csvd_mica_out

# Convert T1w and FLAIR to N4 bias field corrected format using N4BiasFieldCorrection
bids_t1=${bids_dir}/sub-${sub_num}/anat/sub-${sub_num}_T1w.nii.gz
bids_t1_n4=${bids_dir}/sub-${sub_num}/anat/sub-${sub_num}_T1w_n4.nii.gz

bids_falir=${bids_dir}/sub-${sub_num}/anat/sub-${sub_num}_FLAIR.nii.gz
bids_falir_n4=${bids_dir}/sub-${sub_num}/anat/sub-${sub_num}_FLAIR_n4.nii.gz

# N4 bias correction and LST Segmentation
if [ ! -f "${bids_dir}/sub-${sub_num}/wmh/"*lst*.nii.gz ]; then
    echo "Processing subject: ${sub_num}"
    mkdir -p ${bids_dir}/sub-${sub_num}/wmh
    N4BiasFieldCorrection -i $bids_falir -o $bids_falir_n4 -d 3
    N4BiasFieldCorrection -i $bids_t1 -o $bids_t1_n4 -d 3
    # LST
    lst --t1 $bids_t1_n4 --flair $bids_falir_n4 --output ${bids_dir}/sub-${sub_num}/wmh --device cpu >/dev/null
    rm $bids_falir_n4
    rm $bids_t1_n4
else
    echo "Skip subject: ${sub_num}"
fi
