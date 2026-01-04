import ana_utils
import numpy as np
import pandas as pd
import os, sys
from joblib import Parallel, delayed


def par_net_cfc_corr(cm_lab, cm_idx, sub_id, lm_lab):
    """Parallel calculate corr with regional cfc corr"""
    # get sc or cm
    cm_mat = np.load(
        f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm_lab}_{cm_lab}_{ana_utils.ATLAS}.npy"
    )
    # get fc
    fc_mat = np.load(f"{ana_utils.PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy")
    fc_mat = np.power(fc_mat, 2)
    sub_cfc_result = []
    
    # for dynamic cms
    for j in range(cm_mat.shape[0]):
        cm_sub_mat = cm_mat[j]
        np.fill_diagonal(cm_sub_mat, np.nan)
        for yeo_net_name in ana_utils.YEO7_DICT.keys():
            try:
                r_values, p_values = ana_utils.sub_mat_correlatrion(
                    fc_mat[ana_utils.YEO7_DICT[yeo_net_name], :],
                    cm_sub_mat[ana_utils.YEO7_DICT[yeo_net_name], :],
                )
                rsq_values = np.power(r_values, 2)
                rsq_val, p_val = np.mean(rsq_values), np.mean(p_values)

            except Exception as e:
                print(e)
                rsq_val, p_val = np.nan, np.nan

            sub_cfc_result.append(
                {
                    "pid": sub_id,
                    "Group": ana_utils.SUB_ID_DIAG_DICT[sub_id],
                    "cm": cm_lab,
                    "idx": j,
                    "rsq": rsq_val,
                    "pval": p_val,
                    "net": yeo_net_name,
                }
            )
    return sub_cfc_result


def get_cfc_corr(cm_lab, cm_idx, lm_lab):
    """Perform correlation and save csv"""

    # Use subjects from DGN
    sub_id_list = pd.read_csv(f"{ana_utils.PROJ_HOME}/data/sub_list.csv")["pid"].to_list()

    # Main parallel loop
    cfc_results = np.array(
        Parallel(n_jobs=-1)(
            delayed(par_net_cfc_corr)(cm_lab, cm_idx, sub_id, lm_lab)
            for sub_id in sub_id_list
        )
    )

    return cfc_results


# Main
print(f"Run {__file__}", flush=True)
print("-" * 20, flush=True)

# Create the output dir
output_dir = f"{ana_utils.PROJ_HOME}/results/corr_net_lm_cfc"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Delete previous files
os.system(f"find {output_dir}/* -type f -delete")
print(ana_utils.YEO7_DICT, flush=True)
print("-" * 20, flush=True)
cm_arrs = ana_utils.CM_ARR[1:]
merge_tbl = pd.DataFrame()
for cm in cm_arrs:
    cm_lab, cm_idx, _ = cm
    print("Processing: ", cm_lab, flush=True)
    for i, lm_lab in enumerate(ana_utils.LM_LAB):

        cfc_results = get_cfc_corr(cm_lab, cm_idx, lm_lab)

        flattened_list = [item for sublist in cfc_results for item in sublist]
        cfc_results_tbl = pd.json_normalize(flattened_list)

        cfc_results_tbl.to_csv(
            f"{output_dir}/{cm_lab}-cfc_{lm_lab}.csv",
            index=False,
        )
        merge_tbl = pd.concat([merge_tbl, cfc_results_tbl], axis=0)
merge_tbl.to_csv(
            f"{output_dir}/merged-cfc.csv",
            index=False,
        )
import plot_net_cfc