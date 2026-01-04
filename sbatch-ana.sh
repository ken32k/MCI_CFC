#!/bin/bash
###
# @Author: Ken32g ken32k@163.com
# @Date: 2024-06-02 11:59:02
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2025-01-31 19:26:16
# @FilePath: /csvd-sfc/sbatch-ana.sh
# @Description: Complete analystic pipeline
###

source ~/.bashrc

# echo [Info]..... Run svd_proc_mat
# python /public/home/baishw/communication/csvd-sfc/proc_mat/svd_proc_mat.py
# wait

# echo [Info]..... Run svd_proc_fc
# python /public/home/baishw/communication/csvd-sfc/proc_mat/svd_proc_fc.py
# wait

# python /public/home/baishw/communication/csvd-sfc/proc_mat/wm_proc_sc.py

echo [Info]..... Run svd_get_subject_list
rm "${CSVD_DIR}/data/sub_list.csv"
python /public/home/baishw/communication/csvd-sfc/analysis/svd_get_subject_list.py
wait

python /public/home/baishw/communication/csvd-sfc/WMHseg/wm_seg_lst_summary.py
cp "${CSVD_DIR}/data/sub_list_beforePSM_WMHstat.csv" "${CSVD_DIR}/data/sub_list.csv"

# echo [Info]..... Run svd_lm_wmh_vol
# bash /public/home/baishw/communication/csvd-sfc/svd_lm_wmh_vol.sh
# wait

# echo [Info]..... Run plot_sub_info
# python /public/home/baishw/communication/csvd-sfc/analysis/plot_wmh_vol.py
# wait

echo [Info]..... Run corr_glo_lm_cfc
python /public/home/baishw/communication/csvd-sfc/analysis/corr_glo_lm_cfc.py
wait
echo [Info]..... Run lr_glo_lm_cfc
python /public/home/baishw/communication/csvd-sfc/analysis/lr_glo_lm_cfc.py
wait
echo [Info]..... Run plot_dominance
python /public/home/baishw/communication/csvd-sfc/analysis/lr_glo_lm_plot.py
wait

exit
# # python /public/home/baishw/communication/csvd-sfc/analysis/corr_net_lm_cfc.py
# # wait
python /public/home/baishw/communication/csvd-sfc/analysis/lr_net_lm_cfc.py
wait

# python /public/home/baishw/communication/csvd-sfc/analysis/plot_net_cfc.py
# wait

# python /public/home/baishw/communication/csvd-sfc/analysis/lr_net_lm_plot.py
# wait

# # echo [Info]..... Run plot_subject_raw_matrices
# # python /public/home/baishw/communication/csvd-sfc/analysis/plot_subject_raw_matrices.py
# # wait

# echo Done
