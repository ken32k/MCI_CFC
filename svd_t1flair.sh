#!/bin/bash
###
 # @Author: Ken32g ken32k@163.com
 # @Date: 2024-06-29 09:47:29
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-07-01 14:20:32
 # @FilePath: /csvd-sfc/svd_t1flair.sh
 # @Description: T1FLAIR ratio
### 

source $HOME/.bashrc

# mica path
sub_list_path=${CSVD_DIR}/data/sub_list.csv
bids=${CSVD_DIR}/csvd_mica_bids
output_dir=${CSVD_DIR}/csvd_mica_out/t1flair
cd $bids

while IFS= read -r line; do
        array_value=$(echo "$line" | awk -F',' '{print substr($3, 2, length($3)-2)}')
        all_subject_arr+=("$array_value")
done < <(tail -n +2 "$sub_list_path")

# Parallel
echo $all_subject_arr
for subject in ${all_subject_arr[@]}; do
        echo $subject
        (bash /public/home/baishw/apps/KUL_NIS/KUL_T1T2FLAIRMTR_ratio.sh -v -p ${subject: -3} -n 12 &> /dev/null) &
        if (($(wc -w <<<$(jobs -p)) % 4 == 0)); then wait; fi
done

# # Use sbatch to submit subjects
# batch_subjects=("${all_subject_arr[@]:0:24}")
# N=0
# for subject in ${batch_subjects[@]}
# do
#         if [ "$N" -gt 24 ]; then
#                 exit
#         fi
#         N=$((N + 1))
#         echo "[Info] ... Process ${subject: -3}"

#         # sbatch --output="${CSVD_DIR}/log/myelin_${subject}.txt" \
#         #         --error="${CSVD_DIR}/log/myelin_${subject}.err" \
#         #         --job-name="myelin_${subject}" \
#         #         --nodes=1 \
#         #         --ntasks=1 \
#         #         --partition=bme_cpu \
#         #         --cpus-per-task=12 \
#         #         --time=24:00:00 \
#         #         --wrap="bash /public/home/baishw/apps/KUL_NIS/KUL_T1T2FLAIRMTR_ratio.sh -v -p ${subject: -3} -n 12"

# done
