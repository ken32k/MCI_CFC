#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-05-31 22:05:57
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-10-26 16:54:04
# @Description: Submit subject tasks to cluster
###

#!/bin/bash


# export PROJ_HOME=/public/home/baishw/WMH_MCI
export PROJ_HOME=$CSVD_DIR
out_dir=${PROJ_HOME}/csvd_mica_out/micapipe_v0.2.0
bids_dir=${PROJ_HOME}/csvd_mica_bids
max_task_thr=48
for sub_dir in "${out_dir}"/*; do
    (
        sub_id=$(basename "$sub_dir")

        # Extract numeric part starting from the 5th character and remove leading zeros
        sub_num=$(echo "${sub_id:4}" | sed 's/[^0-9]//g' | sed 's/^0*//')
        echo $sub_num
        if [[ "$sub_num" -gt 725 ]]; then
            rsync -av --remove-source-files $sub_dir /public/home/baishw/WMH_MCI/csvd_mica_out/micapipe_v0.2.0
            # rsync -av --remove-source-files $bids_dir /public/home/baishw/WMH_MCI/csvd_mica_bids
            echo $sub_dir
        fi
    ) &
    # Ensure the number of jobs does not exceed the limit (max_task_thr at a time)
    while [ "$(ps aux | grep rsync | wc -l)" -ge $max_task_thr ]; do
        sleep 10
        echo ps aux | grep rsync
        echo "$(date '+%Y-%m-%d %H:%M:%S')"
    done
done

echo "Done."

# rsync -av --remove-source-files /public_bme/data/v-baishw/CSVD /public/home/baishw/WMH_MCI > ~/log/rsync_2041020 &
