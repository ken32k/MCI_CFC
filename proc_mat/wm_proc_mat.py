"""
Author: Ken32g ken32k@163.com
Date: 2024-05-31 23:52:23
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2024-07-02 10:15:33
FilePath: /csvd-sfc/proc_mat/svd_proc_mat.py
Description: Generate communication matrix

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
import wm_matrix as svd_matrix
from scipy.sparse import csr_matrix, diags

from joblib import Parallel, delayed


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
        plt.colorbar().remove()
        plt.axis("off")
    plt.savefig(
        f"{DATA_OUT_DIR}/scs_fig/{sub_id}_scs.png", dpi=300, bbox_inches="tight"
    )
    plt.close()


def reslice_mat(mat):
    """
    Slice the mat, 450 rois to 400 roi
    The output from MICApipe is 450*450 including the subcortical and Cbl
    """
    mat = mat[49:, 49:]
    nrowcol = len(mat) // 2
    mat = np.delete(mat, nrowcol, axis=0)  # delete row at index nrowcol
    mat = np.delete(mat, nrowcol, axis=1)  # delete column at index nrowcol
    return mat


# Insert nan rows and cols if process with err
def remove_all_zero_rows_cols(mat):
    non_zero_idx = np.any(mat != 0.0, axis=0)
    non_zero_mat = mat[non_zero_idx][:, non_zero_idx]
    return non_zero_mat, non_zero_idx


def remove_euclidean_rows_cols(mat, non_zero_idx):
    non_zero_mat = mat[non_zero_idx][:, non_zero_idx]
    return non_zero_mat


def insert_nans(mat, non_zero_idx):
    mat_restore = np.full((NNODE, NNODE), np.nan)
    mat_restore[np.ix_(non_zero_idx, non_zero_idx)] = mat
    return mat_restore


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
    sc_wb = None
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
            sc_reslc = reslice_mat(sc)

            # Threshold
            sc_reslc_thr = svd_matrix.threshold_proportional(sc_reslc, THR_POR)
            sc_reslc_thr = svd_matrix.autofix(sc_reslc_thr)

            # Save the submatrix as a .mat file for subsequent matlab calc
            np.save(
                f"{DATA_OUT_DIR}/scs/{sub_id}_{mat_name}_sc-wei_{ATLAS}.npy",
                sc_reslc_thr,
            )

            scs.append(sc_reslc_thr)

            # If sc-int and sc-wmh does not exist, this allow the following code to be run
            if mat_name == "wb":
                sc_wb = sc_reslc_thr

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
    scs, full_scs = proc_scs(sub_id)

    # If scs lost, then exit
    if not full_scs:
        return

    # Load Euclidean distance matrix
    euc_distance_path = f"{DATA_OUT_DIR}/euc_distance/{sub_id}_euc_{ATLAS}.mat"
    if os.path.exists(euc_distance_path):
        euc_distance_full = sio.loadmat(euc_distance_path)["ed"]
    else:
        print("[Error].... " + sub_id + ": No Euclidean distance exists.")
        return

    for l, sc_wei in enumerate(scs):
        lm_lab = LM_LAB[l]
        # sc_wei = sc.astype(float)
        # sc_bin = networks.binarize_network(sc_wei, 100).astype(int)

        # Create communication matrices and save only with weighted sc

        # Navigation, only wei
        sc_wei, non_zero_idx = remove_all_zero_rows_cols(sc_wei)
        euc_distance = remove_euclidean_rows_cols(euc_distance_full, non_zero_idx)
        try:
            # Need geodesic distance ###### Wed Dec 27 09:46:54 CST 2023
            nav_sr, nav_sr_node, nav_path_len, nav_path_hop, _ = nmtr.navigation_wu(
                euc_distance, sc_wei
            )
            np.save(
                f"{DATA_OUT_DIR}/cm/nav-pl/{sub_id}_{lm_lab}_nav-pl_{ATLAS}.npy",
                [insert_nans((nav_path_len + nav_path_len.T) / 2, non_zero_idx)],
            )
        except Exception as e:
            print(
                f"An error occurred while calculating and saving nav_path_len: {str(e)}"
            )

        # MFPT
        try:
            # remove singularity
            # https://github.com/qsnake/gpaw/blob/d03fb05244a8e3ccc0ca5693e45e030a62e50c4a/gpaw/utilities/tools.py#L121
            # tiny_val = 1e-12
            # sc_wei_add_const = np.where(sc_wei < tiny_val, tiny_val, sc_wei)

            # # Convert mtx to a sparse matrix in CSR format
            # mtx_sparse = csr_matrix(sc_wei)

            # # Add a small constant to the diagonal elements
            # mtx_sparse = mtx_sparse + diags([1e-10] * mtx_sparse.shape[0], 0)

            # # Convert the sparse matrix back to a dense matrix
            # sc_wei_add_const = mtx_sparse.toarray()

            ###### Sat Mar 16 16:12:01 CST 2024
            # https://github.com/Migirditchsv/intuitionMdp/blob/534e719fdce079ba6c21b958013f9f2b772fa397/src/mfpt.py#L79

            mfpt_wei = bc.distance.mean_first_passage_time(sc_wei)

            np.save(
                f"{DATA_OUT_DIR}/cm/mfpt-wei/{sub_id}_{lm_lab}_mfpt-wei_{ATLAS}.npy",
                # nmtr.mean_first_passage_time(sc_wei), # Use nmtr
                [
                    insert_nans((mfpt_wei + mfpt_wei.T) / 2, non_zero_idx)
                ],  # Use bc ###### Mon Mar 11 13:06:40 CST 2024
            )
        except Exception as e:
            print(f"An error occurred while calculating and saving MFPT-wei: {str(e)}")

        # communicability
        try:
            np.save(
                f"{DATA_OUT_DIR}/cm/co-wei/{sub_id}_{lm_lab}_co-wei_{ATLAS}.npy",
                [insert_nans(nmtr.communicability_wei(sc_wei), non_zero_idx)],
            )
        except Exception as e:
            print(
                f"An error occurred while calculating and saving communicability-wei: {str(e)}"
            )

        # SPL, SI and PT
        cm_pls, cm_pts, cm_sis = [], [], []
        gammas = [0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 10]

        for g in gammas:
            # For each gamma, convert SC weight to cost
            sc_cost = np.power(sc_wei, -g)
            sc_cost = np.ascontiguousarray(sc_cost)
            # Shortest path length
            try:
                cm_pls.append(
                    insert_nans(nmtr.distance_wei_floyd(sc_cost)[0], non_zero_idx)
                )
            except:
                cm_pls.append(insert_nans(np.zeros_like(sc_cost), non_zero_idx))
            # Matching index and path transitivity
            # path_transitivity has been modified

            try:
                pt_wei = nmtr.path_transitivity(sc_wei, sc_cost)
                cm_pts.append(insert_nans((pt_wei + pt_wei.T) / 2, non_zero_idx))
            except:
                cm_pts.append(insert_nans(np.zeros_like(sc_cost), non_zero_idx))

            try:
                si_wei = nmtr.search_information(sc_wei, sc_cost)
                cm_sis.append(insert_nans((si_wei + si_wei.T) / 2, non_zero_idx))
            except:
                cm_sis.append(insert_nans(np.zeros_like(sc_cost), non_zero_idx))

            # Save results
            np.save(
                f"{DATA_OUT_DIR}/cm/pl-wei/{sub_id}_{lm_lab}_pl-wei_{ATLAS}.npy",
                cm_pls,
            )

            np.save(
                f"{DATA_OUT_DIR}/cm/pt-wei/{sub_id}_{lm_lab}_pt-wei_{ATLAS}.npy",
                cm_pts,
            )

            np.save(
                f"{DATA_OUT_DIR}/cm/si-wei/{sub_id}_{lm_lab}_si-wei_{ATLAS}.npy",
                cm_sis,
            )

        # Flow graph
        flow_graphs = []
        markovts = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        for mt in markovts:
            try:
                flow_graphs.append(
                    insert_nans(nmtr.flow_graph(sc_wei, t=mt), non_zero_idx)
                )
            except Exception as e:
                flow_graphs.append(insert_nans(np.zeros_like(sc_wei), non_zero_idx))
                print(
                    f"An error occurred while calculating and saving flow graph matrices: {sub_id}, {e}"
                )
        try:
            np.save(
                f"{DATA_OUT_DIR}/cm/fg-wei/{sub_id}_{lm_lab}_fg-wei_{ATLAS}.npy",
                flow_graphs,
            )
        except:
            print(
                f"An error occurred while calculating and saving flow graph matrices: {sub_id}, {e}"
            )
            
def main():
    print(f"Run {__file__}", flush=True)
    global NNODE, ATLAS, TCKNUM, MAT_PREX, THR_POR, DATA_OUT_DIR, MICA_OUT_DIR, LM_LAB

    NNODE = 200  # <<<<<<<<<<<< CHANGE THIS ATLAS, 200 or 400

    ATLAS = "schaefer-" + str(NNODE)
    TCKNUM = "2M"  # <<<<<<<<<<<< CHANGE THIS tcknum
    THR_POR = 0.2  # <<<<<<<<<<<< CHANGE THIS threshold proportion
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
    total_found = len(sub_list)
    total_to_process = len(sub_list)

    print(
        f"Total {total_found} subjects found in the output folder and "
        f"{total_to_process} subjects to process. CSVD project: processing matrices "
        "and deleting existing files."
    )

    # Delete existing files in specified directories
    directories_to_clean = [
        f"{DATA_OUT_DIR}/scs/*",
        f"{DATA_OUT_DIR}/scs-fig/*",
        f"{DATA_OUT_DIR}/cm/*",
        f"{DATA_OUT_DIR}/lms/*",
    ]

    for dir_path in directories_to_clean:
        os.system(f"find {dir_path} -type f -delete")

    # Process matrices in parallel
    Parallel(n_jobs=-1)(delayed(par_proc_mat)(sub_id) for sub_id in sub_list)

if __name__ == "__main__":
    main()
