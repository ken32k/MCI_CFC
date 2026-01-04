#/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-09-16 22:00:28
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-09-17 10:52:46
# @FilePath: /csvd-sfc/parallel_run.sh
# @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
###

export WMH_MCI=/public/home/baishw/WMH_MCI
bids=${WMH_MCI}/csvd_mica_bids
# Parallel
for SUB_DIR in "$bids"/*; do
    (
        # sub_pid=$(basename "$SUB_DIR")
        sub_pid="${SUB_DIR: -3}"
        echo "[Info] $(date +'%Y-%m-%d %H:%M:%S') - Processing $sub_pid"
        bash /public/home/baishw/communication/csvd-sfc/svd_sfc-mica_run.sh $sub_pid 12 2 >$WMH_MCI/logs/${sub_pid}.log
        
        # &>/dev/null
    ) &
    if (($(wc -w <<<$(jobs -p)) % 3 == 0)); then
        echo $(wc -w <<<$(jobs -p))
        echo full wait
        wait
    fi
done
