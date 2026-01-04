
import numpy as np
import pandas as pd
import seaborn as sns
import os
import nibabel as nib


NNODE = 400  # <<<<<<<<<<<< CHANGE THIS ATLAS, 200 or 400

# ATLAS = "schaefer-" + str(NNODE)
TCKNUM = "2M"  # <<<<<<<<<<<< CHANGE THIS tcknum
THR_POR = 0.1  # <<<<<<<<<<<< CHANGE THIS threshold proportion
MAT_PREX = ["wei"]

# Set CSVD paths
HOME = "/public_bme/data/v-baishw/CSVD"
# Matrix output dir
DATA_OUT_DIR = f"{HOME}/data"
sub_data_csv = f"{DATA_OUT_DIR}/sub_info.xlsx"
sub_info = pd.read_excel(sub_data_csv)
SUB_ID_LIST = sub_info.pid.tolist()
print(SUB_ID_LIST)

# Preproc dirs
MICA_OUT_DIR = f"{HOME}/csvd_mica_out/micapipe_v0.2.0"
# MICA_OUT_DIR = "/home/lingbin/ke/svd/raw/bidstest_result/micapipe_v0.2.0"

mica_output_subjects = sorted(os.listdir(MICA_OUT_DIR))

# BIDS dirs
mica_bids_subjects = sorted(os.listdir(f"{HOME}/csvd_mica_bids"))

# Euclidean distance list
euc_dist_subjects = sorted(os.listdir(f"{HOME}/data/euc_distance"))

proc_mat_list = sorted(list(set(SUB_ID_LIST) & set(mica_output_subjects)))
summed_data = []
print(proc_mat_list)
for ATLAS in ["schaefer-200"]:
    for subject in proc_mat_list:
        mica_sub_sc_path = f"{MICA_OUT_DIR}/{subject}/dwi/connectomes-wmh-{str(TCKNUM)}"
        sub_sc_path = f"{mica_sub_sc_path}/{subject}_space-dwi_atlas-{ATLAS}_desc-iFOD2-{str(TCKNUM)}-SIFT2_full-connectome.shape.gii"
        
        if os.path.exists(sub_sc_path):
            try:
                sub_sc = nib.load(sub_sc_path).darrays[0].data
                summed_value = sub_sc.sum(axis=0)  # 按行求和
                summed_data.append([subject+str(ATLAS)] + summed_value.tolist())  # 将每个主体的求和结果与主体名称组合成一个列表，添加到数据列表中
            except Exception as e:
                print(f"An error occurred while processing GIfTI file for subject {subject}: {e}")
        else:
            print(f"GIfTI file for subject {subject} does not exist in the specified path.")

# 创建数据框并保存为CSV
df = pd.DataFrame(summed_data)  # 创建数据框
df.to_csv(f'{HOME}/results/sc_summed_results.csv', index=False)  # 保存为csv文件