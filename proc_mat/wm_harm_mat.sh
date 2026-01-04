###
# @Author: Ken32g ken32k@163.com
# @Date: 2025-01-12 19:13:37
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2025-01-13 21:37:39
# @FilePath: /csvd-sfc/proc_mat/wm_harm_mat.sh
# @Description:
###
export PROJ_HOME=$CSVD_DIR
data_out_dir="${PROJ_HOME}/data"
sub_info_tbl="${data_out_dir}/sub_list_beforePSM_WMHstat.csv"
scs_raw_dir="${data_out_dir}/scs_raw"
scs_out_dir="${data_out_dir}/scs"
covariates="age + sex"

for lm_lab in "wb" "int" "wmh"; do
    echo run $lm_lab
    mat_name_template="${scs_raw_dir}/{s}_${lm_lab}_sc-wei_schaefer-200.npy"
    log_file="${data_out_dir}/harm_log_${lm_lab}.log"
    python /public/home/baishw/communication/csvd-sfc/proc_mat/harmonize_connectomes.py -c $sub_info_tbl -f $mat_name_template -o $scs_out_dir -V "age + sex" -S "Site" --debug --logfile $log_file
done
#  > ${log_file}.txt

# FCs
fcs_raw_dir="${data_out_dir}/fcs_raw"
fcs_out_dir="${data_out_dir}/fcs"
mat_name_template="${fcs_raw_dir}/{s}_${lm_lab}_fc-wei_schaefer-200.npy"
log_file="${data_out_dir}/harm_log_${lm_lab}.log"
python /public/home/baishw/communication/csvd-sfc/proc_mat/harmonize_connectomes.py -c $sub_info_tbl -f $mat_name_template -o $fcs_out_dir -V "age + sex" -S "Site" --debug --logfile $log_file
