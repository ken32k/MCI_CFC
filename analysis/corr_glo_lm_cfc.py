"""
Author: Ken32g ken32k@163.com
Date: 2024-06-01 23:08:10
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2024-07-01 14:38:43
FilePath: /csvd-sfc/analysis/corr_glo_lm_cfc.py
Description: Global univaiate coupling
"""

import ana_utils
import numpy as np
import pandas as pd
import os, sys
from joblib import Parallel, delayed
import statsmodels.api as sm


def par_global_cfc_corr(cm_lab, sub_id, lm_lab):
    """Parallel calculate corr with global cfc corr"""
    # get sc or cm
    cm_mat = np.load(
        f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm_lab}_{cm_lab}_{ana_utils.ATLAS}.npy"
    )
    # get fc
    fc_mat = np.load(
        f"{ana_utils.PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
    )
    sub_cfc_result = []
    # for dynamic cms
    for j in range(cm_mat.shape[0]):
        cm_sub_mat = cm_mat[j]
        try:
            r_values, p_values = ana_utils.sub_mat_correlatrion(fc_mat, cm_sub_mat)
            rsq_values = np.power(r_values, 2)
            rsq_val, p_val = np.mean(rsq_values), np.mean(p_values)

        except Exception as e:
            print(f"error: {e}", flush=True)
            rsq_val, p_val = np.nan, np.nan

        sub_cfc_result.append(
            {
                "pid": sub_id,
                "Group": ana_utils.SUB_ID_DIAG_DICT[sub_id],
                "cm": cm_lab,
                "idx": j,
                "rsq": rsq_val,
                "pval": p_val,
            }
        )
    return sub_cfc_result


def get_cfc_corr(cm_lab, lm_lab):
    """Perform correlation and save csv"""

    # Use subjects from DGN
    sub_id_list = pd.read_csv(
        f"{ana_utils.PROJ_HOME}/data/sub_list.csv", encoding="gbk"
    )["pid"].to_list()

    # Main parallel loop
    cfc_results = np.array(
        Parallel(n_jobs=ana_utils.N_JOBS)(
            delayed(par_global_cfc_corr)(cm_lab, sub_id, lm_lab)
            for sub_id in sub_id_list
        )
    )

    return cfc_results


print(f"Run {__file__}", flush=True)
print("-" * 20, flush=True)

# Create the empty output dir
output_dir = f"{ana_utils.PROJ_HOME}/results/corr_glob_lm_cfc"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.system(f"find {output_dir}/* -type f -delete")


cm_arr = ana_utils.CM_ARR[1:]

for cm in cm_arr:
    cm_lab = cm[0]
    for lm_lab in ana_utils.LM_LAB:
        print(f"processing {cm_lab}", flush=True)
        cfc_results = get_cfc_corr(cm_lab, lm_lab)

        flattened_list = [item for sub_id_list in cfc_results for item in sub_id_list]
        cfc_results_tbl = pd.json_normalize(flattened_list)

        cfc_results_tbl.to_csv(
            f"{output_dir}/{cm_lab}-cfc_{lm_lab}.csv",
            index=False,
        )

import plot_all_cms
import plot_dyn_cms
