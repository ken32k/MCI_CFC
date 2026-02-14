#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-08-28 13:00:31
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-10-03 14:35:14
# @FilePath: //sbatch_wmh_seg.sh
# @Description: Sbatch WMH segmentation with LST or synthseg
###

# load fs_dev
export PROJ_HOME=$MCI_CFC_DIR
export FREESURFER_HOME=/public_bme/data/v-baishw/app/freesurfer_dev
source $FREESURFER_HOME/SetUpFreeSurfer.sh

export FSF_OUTPUT_FORMAT=nii
# export SUBJECTS_DIR=${PROJ_HOME}/mica_out/fastsurfer
export MNI_DIR=${FREESURFER_HOME}/mni

# mica path
bids=${PROJ_HOME}/mica_bids
out=${PROJ_HOME}/mica_out
tmp=${PROJ_HOME}/mica_temp

squeue_prex=w
for sub_dir in ${bids}/*; do
    (
        # target_file=${sub_dir}/anat/*WMH_seg.nii.gz
        target_file="${bids}/${sub_id}/wmh/"*lst*.nii.gz
        sub_id=$(basename ${sub_dir})
        if [ ! -f $target_file ]; then
            sbatch --output="${sub_dir}/sbatch.txt" \
                --error="${sub_dir}/sbatch.err" \
                --job-name="${squeue_prex}$(basename $sub_dir)" \
                --nodes=1 \
                --ntasks=1 \
                --partition=bme_cpu \
                --cpus-per-task=10 \
                --time=00:50:00 \
                /public/home/baishw/communication//WMHseg/mci_wmh_seg_lst.sh ${sub_id:4}

            echo "Submit $sub_dir"
        else
            echo "Skip $sub_dir"
        fi
    ) &
    # Check if the squeue greater than threshold
    while [ $(squeue -u baishw | grep ${squeue_prex}sub | wc -l) -ge 12 ]; do
        sleep 2
        echo $(date '+%Y-%m-%d %H:%M:%S')
    done
done
echo done
