import ana_utils
import numpy as np
import pandas as pd
import seaborn as sns
import os, sys
import scipy as sp
import nilearn as nil
from nilearn import plotting
from scipy.stats import f_oneway
import statsmodels.api as sm
import matplotlib.pyplot as plt
from joblib import Parallel, delayed


def par_global_cfc_lrm(cm_list, sub_id, lm_lab):
    """Parallel calculate corr with global lrm"""

    # Get communication model matrices
    cm_mats, dominance = [], []
    for cm in cm_list:
        cm_lab, cm_idx = cm[0], cm[1]
        cm_mat = np.load(
            f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm_lab}_{cm_lab}_{ana_utils.ATLAS}.npy"
        )
        cm = cm_mat[cm_idx][ana_utils.TRIU_IDX]
        cm_mats.append(cm)
    # Get fc
    fc_mat = np.load(
        f"{ana_utils.PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
    )[ana_utils.TRIU_IDX]
    # Get rsq, pval and dominance
    try:
        rsq_val, p_val, dominance, rsq_val_domin, lenan = ana_utils.sub_mat_linear_reg(
            fc_mat, np.array(cm_mats)
        )
    except Exception as e:
        rsq_val, p_val, lenan = np.nan, np.nan, np.nan
        dominance = [np.nan for i in range(5)]

    sub_lrm_result = {
        "pid": sub_id,
        "Group": ana_utils.SUB_ID_DIAG_DICT[sub_id],
        "cm": "lrm",
        "idx": 0,
        "rsq": rsq_val_domin,
        "rsq_lrm": rsq_val,
        "pval": p_val,
        # "lenan": lenan,
    }
    # Add dominance
    for i, cm in enumerate(cm_list):
        cm_lab, cm_idx = cm[0], cm[1]
        sub_lrm_result[cm_lab] = dominance[i]
    return sub_lrm_result


def get_cfc_lrm(cm_list, lm_lab):
    """Perform correlation and save csv"""
    # Get subject list
    sub_cm_list = pd.read_csv(f"{ana_utils.PROJ_HOME}/data/sub_list.csv")[
        "pid"
    ].to_list()

    # Main parallel loop
    # cfc_results = np.array(
    #     Parallel(n_jobs=6)(
    #         delayed(par_global_cfc_lrm)(cm_list, sub_id, lm_lab)
    #         for sub_id in sub_cm_list
    #     )
    # )

    # Main loop
    cfc_results = np.array([])
    for idx, sub_id in enumerate(sub_cm_list):
        print(lm_lab, idx, sub_id,flush=True)
        cfc_results = np.append(
            cfc_results, par_global_cfc_lrm(cm_list, sub_id, lm_lab)
        )
    return cfc_results


print(f"Run {__file__}", flush=True)
print("-" * 20, flush=True)

# Create the output dir
output_dir = f"{ana_utils.PROJ_HOME}/results/lr_glob_lm_cfc"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# # Delete existing files
os.system(f"find {output_dir}/* -type f -delete")

cm_list = ana_utils.CM_ARR[1:]

# Get RSQ with Linear regression model
for i, lm_lab in enumerate(ana_utils.LM_LAB):
    print(lm_lab, flush=True)
    cfc_results = get_cfc_lrm(cm_list, lm_lab)
    # print(cfc_results, flush=True)
    flattened_list = [item for item in cfc_results]
    cfc_results_tbl = pd.json_normalize(flattened_list)
    cfc_results_tbl.to_csv(
        f"{output_dir}/lrcms-cfc_{lm_lab}.csv",
        index=False,
    )

    for cm in cm_list:
        cm_lab, cm_idx = cm[0], cm[1]
        cfc_results_tbl.loc[:, "rsq"] = cfc_results_tbl[cm_lab]
        cm_dom_tbl = cfc_results_tbl[["pid", "Group", "cm", "idx", "rsq"]]
        cm_dom_tbl.loc[:, "idx"] = cm_idx
        cm_dom_tbl.to_csv(
            f"{output_dir}/{cm_lab}-cfc_{lm_lab}.csv",
            index=False,
        )
import lr_glo_lm_plot
