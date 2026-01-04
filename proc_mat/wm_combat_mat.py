"""
Author: Ken32g ken32k@163.com
Date: 2025-01-05 21:03:10
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2025-01-05 21:06:18
FilePath: /csvd-sfc/proc_mat/wm_combat_mat.py
Description: 
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas as pd
from neuroCombat import neuroCombat


def get_upper_log_mat(raw_mat):
    upper_tri_indices = np.triu_indices(raw_mat.shape[0], k=1)
    upper_tri_values = raw_mat[upper_tri_indices]
    upper_tri_values[upper_tri_values == 0] = 1e-12
    upper_tri_values_log = np.log(upper_tri_values)
    return upper_tri_indices, upper_tri_values_log


def main():
    print(f"Run {__file__}", flush=True)
    print("=======Processing Combat=========")

    nnode = 200  # <<<<<<<<<<<< CHANGE THIS ATLAS, 200 or 400
    thr_por = 0.25  # <<<<<<<<<<<< CHANGE THIS threshold proportion
    lm_labs = ["wb", "wmh", "int"]
    # Set CSVD paths
    HOME = "/public/home/baishw/WMH_MCI"

    # Matrix output dir
    DATA_OUT_DIR = f"{HOME}/data"
    sub_info_tbl = pd.read_csv(f"{DATA_OUT_DIR}/sub_list_beforePSM_WMHstat.csv")
    sub_info_tbl = sub_info_tbl[
        sub_info_tbl.pid.isin(
            set([item[:7] for item in os.listdir(f"{DATA_OUT_DIR}/scs_raw/")])
        )
    ]
    sub_list = sub_info_tbl.pid.to_list()

    # print(sub_list, flush=True)
    center_lab = []

    # Iterate over each subject ID
    center_lab = [1 if int(sub_id.split("-")[1]) <= 630 else 2 for sub_id in sub_list]
    # Display the resulting list
    print(center_lab, flush=True)

    # Process SC
    for lm_lab in lm_labs:
        print(f"------Processing {lm_lab}------", flush=True)
        multicenter_data = []
        for sub_id in sub_list:
            sub_mat_path = (
                f"{DATA_OUT_DIR}/scs_raw/{sub_id}_{lm_lab}_sc-wei_schaefer-200.npy"
            )
            # Perform log transformation
            raw_sc_mat = np.load(sub_mat_path)
            print(sub_id, np.max(raw_sc_mat), flush=True)
            upper_tri_indices, upper_tri_values_log = get_upper_log_mat(raw_sc_mat)
            multicenter_data.append(upper_tri_values_log)

        multicenter_tbl = pd.DataFrame(np.column_stack(multicenter_data))

        covars = {
            "batch": center_lab,
            "age": sub_info_tbl.age.to_list(),
            "sex": sub_info_tbl.sex.to_list(),
            "logwmh": sub_info_tbl.logwmh.to_list(),
        }
        covars = pd.DataFrame(covars)

        categorical_cols = ["sex"]
        batch_col = "batch"
        # print(multicenter_tbl, flush=True)
        multicenter_tbl.to_csv(f"{DATA_OUT_DIR}/combat/sc_pre_{lm_lab}.csv")
        # Harmonization step:
        data_combat = neuroCombat(
            dat=multicenter_tbl,
            covars=covars,
            batch_col=batch_col,
            categorical_cols=categorical_cols,
            ref_batch="1",
            mean_only=True,
        )["data"]
        pd.DataFrame(data_combat).to_csv(f"{DATA_OUT_DIR}/combat/sc_data_combat_{lm_lab}.csv",header=False)
        # Write harmonized matrix for each subject
        for idx, sub_id in enumerate(sub_list):

            full_mat = np.zeros((nnode, nnode))
            upper_tri_indices = np.triu_indices(nnode, k=1)

            # Assign the sliced data_combat values to the upper triangular part
            full_mat[upper_tri_indices] = np.exp(data_combat[:, idx])

            # Transpose to flip upper to lower triangular
            sqr_mat = full_mat + full_mat.T

            raw_sc_mat = np.load(
                f"{DATA_OUT_DIR}/scs_raw/{sub_id}_{lm_lab}_sc-wei_schaefer-{nnode}.npy"
            )
            raw_sc_mat[raw_sc_mat > 0] = 1
            sqr_mat = sqr_mat * raw_sc_mat

            # Threshold SC
            sqr_mat_thr = wm_matrix.threshold_proportional(sqr_mat, thr_por)

            # np.fill_diagonal(sqr_mat, 0)
            np.save(
                f"{DATA_OUT_DIR}/scs/{sub_id}_{lm_lab}_sc-wei_schaefer-{nnode}.npy",
                sqr_mat_thr,
            )
            print(
                f"{sub_id} | Negative values: {np.any(sqr_mat <0)}, None-zero: {np.count_nonzero(sqr_mat_thr) / (nnode**2) * 100} % from {np.count_nonzero(sqr_mat) / (nnode**2) * 100} %",
                flush=True,
            )

    # FC
    print(f"------Processing FC------", flush=True)
    multicenter_data = []
    for sub_id in sub_list:
        print()
        sub_mat_path = f"{DATA_OUT_DIR}/fcs_raw/{sub_id}_fc-wb_schaefer-{nnode}.npy"
        mat = np.load(sub_mat_path)
        upper_tri_indices = np.triu_indices(mat.shape[0], k=1)
        upper_tri_values = mat[upper_tri_indices]
        multicenter_data.append(upper_tri_values)

    multicenter_tbl = pd.DataFrame(np.column_stack(multicenter_data))

    covars = {
        "batch": center_lab,
        "age": sub_info_tbl.age.to_list(),
        "sex": sub_info_tbl.sex.to_list(),
        "logwmh": sub_info_tbl.logwmh.to_list(),
    }
    covars = pd.DataFrame(covars)

    categorical_cols = ["sex"]

    batch_col = "batch"
    print(multicenter_tbl, flush=True)
    multicenter_tbl.to_csv(f"{DATA_OUT_DIR}/combat/fc_pre.csv")
    # Harmonization step:
    data_combat = neuroCombat(
        dat=multicenter_tbl,
        covars=covars,
        batch_col=batch_col,
        categorical_cols=categorical_cols,
        ref_batch="1",
        mean_only=False,
    )["data"]
    pd.DataFrame(data_combat).to_csv(
        f"{DATA_OUT_DIR}/combat/fc_data_combat.csv", header=False
    )

    # Write harmonized matrix for each subject
    for idx, sub_id in enumerate(sub_list):

        full_matrix = np.zeros((nnode, nnode))
        upper_tri_indices = np.triu_indices(nnode, k=1)

        # Assign the sliced data_combat values to the upper triangular part
        full_matrix[upper_tri_indices] = data_combat[:, idx]

        # Transpose to flip upper to lower triangular
        lower_tri_matrix = full_matrix + full_matrix.T
        np.fill_diagonal(lower_tri_matrix, 0)
        np.save(
            f"{DATA_OUT_DIR}/fcs/{sub_id}_fc-wb_schaefer-{nnode}.npy",
            lower_tri_matrix,
        )


if __name__ == "__main__":
    main()
