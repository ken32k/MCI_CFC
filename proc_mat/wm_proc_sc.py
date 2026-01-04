"""
Author: Ken32g ken32k@163.com
Date: 2024-05-31 23:52:23
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2024-07-02 10:15:33
FilePath: /csvd-sfc/proc_mat/svd_proc_mat.py
Description: Process structural connectomes

 Use miniconda 3.11.5 base on BME cluster
 For jupyter: sometimes the vscode has a delay, just wait for 3.11.5 shown
"""

import os, datetime, pdb, pickle
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas as pd

from tqdm import tqdm
from nilearn import plotting
from netneurotools import modularity, networks
from netneurotools import metrics as nmtr
import brainconn as bc
import wm_matrix
from scipy.sparse import csr_matrix, diags

from joblib import Parallel, delayed
import argparse


def flip_upper_triangle(matrix):
    # new_matrix = np.copy(matrix)
    # upper_triangle_indices = np.triu_indices(matrix.shape[0], 1)
    # new_matrix[(upper_triangle_indices[1], upper_triangle_indices[0])] = matrix[
    #     upper_triangle_indices
    # ]

    flipped_mat = np.transpose(matrix)

    # Add the flipped matrix to the original matrix
    result = matrix + flipped_mat

    # Divide the matrix element-wise by 2
    new_matrix = result / 2.0
    return new_matrix


def save_mat_image(matrices, sub_id, cm="Purples"):
    """
    Save a png picture of the SCs
    Binarize weighted matrix
    """
    matrices_bin = [(m > 0).astype(int) for m in matrices]
    matrices = matrices + matrices_bin
    # Visualize the data
    for i in range(6):
        plt.subplot(2, 3, i + 1)
        cm = "OrRd_r" if np.any(np.sum(matrices[i]) == 0) else "Purples"
        plt.imshow(matrices[i], cmap=cm)
        # plt.colorbar().remove()
        plt.axis("off")
    plt.savefig(
        f"{DATA_OUT_DIR}/scs_fig/{sub_id}_scs.png", dpi=300, bbox_inches="tight"
    )
    plt.close()


def reslice_mica_mat(sc_mat):
    """
    Slice the sc_mat, 450 rois to 400 roi
    The output from MICApipe is 450*450 including the subcortical and Cbl
    """
    sc_mat = sc_mat[49:, 49:]
    nrowcol = len(sc_mat) // 2
    sc_mat = np.delete(sc_mat, nrowcol, axis=0)  # delete row at index nrowcol
    sc_mat = np.delete(sc_mat, nrowcol, axis=1)  # delete column at index nrowcol
    return sc_mat


# Process SCs
def proc_scs(sub_id):
    """
    Process structural connectivity.
    return the whole brain SC

    Args:
        sub_id (string): subject ID
    """
    full_scs = 1
    scs = []
    sc_exist = True
    for mat_name in LM_LAB:
        # sub_id sc
        mica_sub_sc_path = (
            f"{MICA_OUT_DIR}/{sub_id}/dwi/connectomes-{mat_name}-{str(TCKNUM)}"
        )

        # Connectome path
        sub_sc = f"{mica_sub_sc_path}/{sub_id}_space-dwi_atlas-{ATLAS}_desc-iFOD2-{str(TCKNUM)}-SIFT2_full-connectome.shape.gii"

        # Gii2mat
        if os.path.exists(sub_sc):
            # Load the GIFTI file
            sc = nib.load(sub_sc).darrays[0].data
            sc = flip_upper_triangle(sc)
            sc_reslc = reslice_mica_mat(sc)

            # Threshold
            # sc_reslc_thr = wm_matrix.threshold_proportional(sc_reslc, THR_POR)
            # sc_reslc_thr = wm_matrix.autofix(sc_reslc_thr)

            ###### Fri Jan 31 11:29:39 CST 2025
            # threshold after combat
            sc_reslc_thr = sc_reslc
            if np.max(sc_reslc_thr) < 100:
                # Save the submatrix as a .mat file for subsequent matlab calc
                np.save(
                    f"{DATA_OUT_DIR}/scs_raw/{sub_id}_{mat_name}_sc-wei_{ATLAS}.npy",
                    sc_reslc_thr,
                )

                scs.append(sc_reslc_thr)
            else:
                sc_exist = False
                scs.append(np.zeros((NNODE, NNODE)))
                full_scs = 0

        else:
            # print("[Error].... " + sub_id + ":no dir" + mica_sub_sc_path, flush=True)
            sc_exist = False
            scs.append(np.zeros((NNODE, NNODE)))
            full_scs = 0
            pass

    if not sc_exist:
        print("[Error].... " + sub_id + ":SC file not exist", flush=True)

    # Save the figure for wb, wmh and int SCs
    save_mat_image(scs, sub_id)
    return scs, full_scs


def par_proc_mat(sub_id):
    print("[Info] .... Proc: ", sub_id, flush=True)

    # Process scs
    _, full_scs = proc_scs(sub_id)

    # If any of scs lost, then exit
    if not full_scs:
        return


def main():
    print(f"Run {__file__}", flush=True)
    global NNODE, ATLAS, TCKNUM, MAT_PREX, THR_POR, DATA_OUT_DIR, MICA_OUT_DIR, LM_LAB

    # Parameters
    NNODE = 200  # <<<<<<<<<<<< CHANGE THIS ATLAS, 200 or 400

    ATLAS = "schaefer-" + str(NNODE)
    TCKNUM = "2M"  # <<<<<<<<<<<< CHANGE THIS tcknum
    THR_POR = 0.3  # <<<<<<<<<<<< CHANGE THIS threshold proportion
    MAT_PREX = ["wei"]
    LM_LAB = ["wb", "wmh", "int"]

    # Set CSVD paths
    HOME = "/public/home/baishw/WMH_MCI"
    # Directory of the micapipe output files
    MICA_OUT_DIR = f"{HOME}/csvd_mica_out/micapipe_v0.2.0"
    # Matrix output dir
    DATA_OUT_DIR = f"{HOME}/data"

    sub_list = sorted(os.listdir(MICA_OUT_DIR))

    # Summary of subject IDs
    total_to_process = len(sub_list)

    print(
        f"Totally {total_to_process} subjects to process. "
        "CSVD project: processing matrices and deleting existing files."
        "=======Delete existing files======="
    )

    # Delete existing files in specified directories
    directories_to_clean = [
        f"{DATA_OUT_DIR}/scs_raw/*",  # directory for raw SCs
        f"{DATA_OUT_DIR}/scs-fig/*",  # directory for raw SCs figures
        # f"{DATA_OUT_DIR}/lms/*",
    ]

    for dir_path in directories_to_clean:
        os.system(f"find {dir_path} -type f -delete || true")
    print("=======Processing SC=========")
    # Process matrices in parallel
    Parallel(n_jobs=args.n_jobs)(delayed(par_proc_mat)(sub_id) for sub_id in sub_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process matrices in parallel.")

    parser.add_argument(
        "--n-jobs",
        type=int,
        default=12,
        help="Number of parallel jobs to run (default: 12)",
    )

    args = parser.parse_args()
    main()
