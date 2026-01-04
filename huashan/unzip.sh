#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-08-28 11:04:16
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-10-02 19:43:45
# @FilePath: /csvd-sfc/huashan/svd_wmh-seg.sh
# @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
###

cd /public_bme/share/HUASHAN/PET_center/MRI/
for zipfile in AFM*.zip; do
    (
        if [[ -f "$zipfile" ]]; then       # 确保是文件
            dirname="${zipfile%.zip}"      # 获取去掉扩展名的文件名
            mkdir -p "$dirname"            # 创建目录
            unzip "$zipfile" -d "$dirname" # 解压到该目录
        fi
    ) &
    if (($(wc -w <<<$(jobs -p)) % 24 == 0)); then
        echo $(wc -w <<<$(jobs -p))
        echo full wait
        wait
    fi
done
