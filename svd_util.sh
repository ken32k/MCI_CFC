###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-11-26 11:19:33
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-12-11 21:44:57
# @FilePath: /csvd-sfc/svd_util.sh
# @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
###
# countt=0; for sub in ./*;do count=$(find $sub/dwi -name connect* -type d | wc -l); if [[ "$count" -eq 1 ]]; then echo $sub:${count}; let countt=$countt+1; fi; done; echo $countt

# countt=0; for sub in ./*;do count=$(find $sub/func -name *schaefer-200_desc-FC.shape.gii -type d | wc -l); if [[ "$count" -eq 1 ]]; then echo $sub:${count}; let countt=$countt+1; fi; done; echo $countt

#!/bin/bash

# 遍历当前目录下的所有子目录
for sub_dir in */; do
    # 进入子目录
    sub_id=$(basename "$sub_dir")
    # Extract numeric part starting from the 5th character and remove leading zeros
    sub_num=$(echo "${sub_id:4}" | sed 's/[^0-9]//g' | sed 's/^0*//')
    if [ "$sub_num" -gt 630 ]; then
        # echo $sub_num
        chckfile=${sub_dir}/dwi/connectomes-wb-2M/${sub_id}_space-dwi_atlas-schaefer-200_desc-iFOD2-2M-SIFT2_full-connectome.shape.gii

        if [ -f "$chckfile" ]; then
            # 获取文件的最后访问时间
            access_time=$(stat -c %X "$chckfile")

            # 将 2024-12-11 转换为 Unix 时间戳
            cutoff_time=$(date -d "2024-12-11" +%s)

            if [ "$access_time" -lt "$cutoff_time" ]; then
                echo "Deleting folder: ${sub_dir}/dwi/connectomes-wb-2M"
                rm -rf "${sub_dir}/dwi/connectomes-wb-2M"
            else
                echo "Access time is after the cutoff date, skipping deletion."
            fi
        else
            echo "File does not exist: $chckfile"
        fi
    fi
done
