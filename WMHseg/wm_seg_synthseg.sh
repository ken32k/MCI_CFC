#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-08-27 19:24:11
# @LastEditors: Ken32g ken32k@163.com
# @LastEditTime: 2024-09-26 22:13:06
# @FilePath: /csvd-sfc/svd_wmh_seg.sh
# @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
###
#
sub_dir=$1
nthr=$2
export PROJ_HOME=$CSVD_DIR
# load fs_dev
export FREESURFER_HOME=/public_bme/data/v-baishw/app/freesurfer_dev
source $FREESURFER_HOME/SetUpFreeSurfer.sh

export FSF_OUTPUT_FORMAT=nii
# export SUBJECTS_DIR=${PROJ_HOME}/csvd_mica_out/fastsurfer
export MNI_DIR=${FREESURFER_HOME}/mni

# mica path
bids_dir=${PROJ_HOME}/csvd_mica_bids

cd ${sub_dir}/anat
sub_id=$(basename ${sub_dir})
echo ${sub_id}
flair_path=$sub_dir/anat/${sub_id}_FLAIR.nii.gz
wmh_out_path=$sub_dir/anat/${sub_id}_WMH_seg.nii.gz

# check if the flair file exists and the wmh segmentation file doesn't exist
if
    [ ! -f $wmh_out_path ] &
    [ -f $flair_path ]
then
    echo "Start to process $sub_id"
    mri_WMHsynthseg --i $flair_path --o $wmh_out_path --threads $nthr --save_lesion_probabilities
else
    echo "$sub_id's WMH already exists, skip it."
fi
