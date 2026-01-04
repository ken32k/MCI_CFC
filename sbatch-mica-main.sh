#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-05-31 22:05:57
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-12-15 09:03:06
# @Description: Submit subject tasks to cluster
###
max_task_thr=20
cpu_per_task=6
# export PROJ_HOME=/public/home/baishw/WMH_MCI
export PROJ_HOME=$CSVD_DIR
bids_dir=${PROJ_HOME}/csvd_mica_bids
out_dir=${PROJ_HOME}/csvd_mica_out

# Alternatively
# Get the subject list from the PSM csv file
# sub_list_tbl=${PROJ_HOME}/data/sub_list_psm2515.csv

# # Get the column "pid" (3rd column) in the CSV into an array, ignoring the header
# sub_list=($(awk -F',' 'NR>1 {print $3}' "$sub_list_tbl"))
# echo "${sub_list[@]}"

for sub_dir in "${bids_dir}"/*; do
    (
        sub_id=$(basename "$sub_dir")
        # Extract numeric part starting from the 5th character and remove leading zeros
        sub_num=$(echo "${sub_id:4}" | sed 's/[^0-9]//g' | sed 's/^0*//')

        # Skip if sub_num is less than 630
        if [[ "$sub_num" -lt 630 ]]; then
            continue
        fi

        # Determine MICA version based on sub_num
        mica_version=1
        if [[ "$sub_num" -gt 630 ]]; then
            mica_version=2
        fi
        count_files=0

        for i in {1..3}; do
            target_file=${out_dir}/micapipe_v0.2.0/${sub_id}/dwi/connectomes-$(case $i in 1) echo "wb" ;; 2) echo "int" ;; 3) echo "wmh" ;; esac)-2M/${sub_id}_space-dwi_atlas-schaefer-200_desc-iFOD2-2M-SIFT2_full-connectome.shape.gii
            if [ -f "$target_file" ]; then
                ((count_files++))
            fi
        done
        echo ${sub_id} ${count_files}

        # Main check: File existence and sub_id presence in the list
        if [[ "$count_files" -lt 3 ]]; then
            rm -rf ${out_dir}/micapipe_v0.2.0/${sub_id}/dwi/connectomes-wb-2M
            rm -rf ${out_dir}/micapipe_v0.2.0/${sub_id}/dwi/connectomes-int-2M
            rm -rf ${out_dir}/micapipe_v0.2.0/${sub_id}/dwi/connectomes-wmh-2M
            
            sbatch --output="${sub_dir}/sbatch-mica.txt" \
                --error="${sub_dir}/sbatch-mica.err" \
                --job-name="m$(basename "$sub_dir")" \
                --nodes=1 \
                --ntasks=1 \
                --partition=bme_cpu \
                --cpus-per-task=${cpu_per_task} \
                --time=12:00:00 \
                /public/home/baishw/communication/csvd-sfc/svd_sfc-mica_run.sh "${sub_num}" ${cpu_per_task} ${mica_version}

            echo "Submitted $sub_dir at $(date '+%Y-%m-%d %H:%M:%S') $count_files"
        else
            echo "Skipped $sub_dir"
            rm -rf ${out_dir}/micapipe_v0.2.0/${sub_id}/dwi/eddy
        fi

    ) &

    # Ensure the number of jobs does not exceed the limit (max_task_thr at a time)
    while [ "$(squeue -u baishw | grep msub | wc -l)" -ge $max_task_thr ]; do
        sleep 10
        # echo "$(date '+%Y-%m-%d %H:%M:%S')"
    done
done

echo "Done."
