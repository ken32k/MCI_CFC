import ana_utils
import numpy as np
import pandas as pd
import os
from joblib import Parallel, delayed


def check_matrix_nan(matrix):
    """
    Check all-zeros, nan and inf for the matrix
    Args:
        matrix (numpy.ndarray): matrix
    Returns:
         (Boollean): contain or not
    """

    has_nan_or_inf = np.isnan(matrix).any() or np.isinf(matrix).any()
    all_zeros = np.all(matrix == 0)

    if has_nan_or_inf or all_zeros:
        return True
    else:
        return False


def load_proc_cm(cm_idx, sub_id, cm, lm_lab):
    cm_lab = cm[0]
    fc_mat_path = f"{ana_utils.PROJ_HOME}/data/fcs_raw/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
    cm_mat_path = f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm_lab}_{cm_lab}_{ana_utils.ATLAS}.npy"
    if os.path.exists(cm_mat_path) and os.path.exists(fc_mat_path):
        cm_mat = np.load(cm_mat_path)[cm[1]]
        # Set the diagonals to 0
        np.fill_diagonal(cm_mat, 0)

        # Replace the inf to 0 in nav_pl
        cm_mat = np.where(np.isinf(cm_mat), 0, cm_mat)

        # Check and remove all zero and nan-containing cms
        if not check_matrix_nan(cm_mat):
            return cm_mat
        else:
            return []
    else:
        return []


# Main
print(f"Run {__file__}", flush=True)
print(f"================", flush=True)

print(len(ana_utils.SUB_ID_LIST))
out_dir = f"{ana_utils.PROJ_HOME}/data"
sub_list = ana_utils.SUB_ID_LIST
print(sub_list, flush=True)
print(out_dir, flush=True)
# Count subjects for each stage
MICA_OUT_DIR = f"{ana_utils.PROJ_HOME}/csvd_mica_out/micapipe_v0.2.0"
BIDS_DIR = f"{ana_utils.PROJ_HOME}/csvd_mica_bids"
bids_count = 0
mica_count = 0
bids_list = []
mica_list = []

# Remove existing sub_list file
if os.path.exists(f"{ana_utils.PROJ_HOME}/data/sub_list.csv"):
    os.remove(f"{ana_utils.PROJ_HOME}/data/sub_list.csv")
# Count subjects with both DWI and BOLD files
for sub in sub_list:
    if (
        # os.path.exists(f"{BIDS_DIR}/{sub}/func/{sub}_dir-AP_rest_bold.nii.gz")
        # and os.path.exists(f"{BIDS_DIR}/{sub}/dwi/{sub}_dir-AP_dwi.nii.gz")
        # and os.path.exists(f"{BIDS_DIR}/{sub}/anat/{sub}_WMH.nii.gz")
        os.path.exists(f"{BIDS_DIR}/{sub}/anat/{sub}_FLAIR.nii.gz")
    ):
        bids_count += 1
        bids_list.append(sub)


# Count subjects with both SC and BOLD FC
# for sub in sub_list:
#     if os.path.exists(
#         f"{MICA_OUT_DIR}/{sub}/dwi/connectomes-wmh-2M/{sub}_surf-fsLR-5k_desc-iFOD2-2M-SIFT2_full-edgeLengths.shape.gii"
#     ) and os.path.exists(
#         f"{MICA_OUT_DIR}/{sub}/func/desc-se_dir-AP_rest_bold/surf/{sub}_surf-fsLR-32k_atlas-schaefer-200_desc-FC.shape.gii"
#     ):
#         mica_count += 1
#         mica_list.append(sub)

# Count CMs
for sub in sub_list:
    if os.path.exists(
        f"{MICA_OUT_DIR}/{sub}/dwi/connectomes-wmh-2M/{sub}_space-dwi_atlas-schaefer-200_desc-iFOD2-2M-SIFT2_full-connectome.shape.gii"
    ) and os.path.exists(
        f"{MICA_OUT_DIR}/{sub}/func/desc-se_dir-AP_rest_bold/surf/{sub}_atlas-schaefer-200_desc-FC.shape.gii.shape.gii"
    ):
        mica_count += 1
        mica_list.append(sub)

print(
    "Folders in BIDS:",
    bids_count,
    "Folders in MICA:",
    mica_count,
    "Folders in FILT:",
    len(pd.read_csv(f"{out_dir}/filt_sub_info.csv", encoding="gbk")),
)


# Filter those who have both wb and int CMs
sub_id_filt = []
x_cm = np.empty((0, 200, 200, len(ana_utils.CM_ARR[1:])))
y_fc = np.empty((0, 200, 200))
for sub_id in sub_list:
    lack_mat = []
    isfull = True
    cm_mats = np.zeros((200, 200, len(ana_utils.CM_ARR[1:])))
    for cm_idx, cm in enumerate(ana_utils.CM_ARR[1:]):
        cm_mat = load_proc_cm(cm_idx, sub_id, cm, "wb")
        if len(cm_mat):
            cm_mats[:, :, cm_idx] = cm_mat
        else:
            isfull = False
            lack_mat = [*lack_mat, cm[0]]
            # print(sub_id, cm[0], flush=True)
            # break

    if isfull:
        fc_mat = np.load(
            f"{ana_utils.PROJ_HOME}/data/fcs_raw/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
        )
        # np.fill_diagonal(fc_mat, 0)
        # x_cm = np.append(x_cm, cm_mats[np.newaxis], axis=0)
        # y_fc = np.append(y_fc, fc_mat[np.newaxis], axis=0)
        if (
            ana_utils.SUB_ID_AGE_DICT[sub_id] >= 50
            or ana_utils.SUB_ID_AGE_DICT[sub_id] <= 80
        ):
            sub_id_filt.append(sub_id)
            print("-" * 12, sub_id, ": In", flush=True)
        else:
            print("-" * 12, sub_id, ": Ex--age", flush=True)
    else:
        print("-" * 12, sub_id, ": Ex--Lack", lack_mat, flush=True)


filtered_sub_id_diag_dict = {
    key: ana_utils.SUB_ID_DIAG_DICT[key] for key in sub_id_filt
}

value_counts = {
    value: sum(1 for v in filtered_sub_id_diag_dict.values() if v == value)
    for value in set(filtered_sub_id_diag_dict.values())
}

print(
    f"anautil subject list:",
    len(ana_utils.SUB_ID_LIST),
    "Final subjects: ",
    value_counts,
    "final len ",
    str(len(x_cm)),
    flush=True,
)

exclude = ["sub-304", "sub-350", "sub-289", "sub-045"]
sub_id_filt = [sub for sub in sub_id_filt if sub not in exclude]

filt_sub_tbl = pd.DataFrame({"pid": sub_id_filt})
filt_sub_tbl["Group"] = filt_sub_tbl["pid"].map(ana_utils.SUB_ID_DIAG_DICT)
filt_sub_tbl["age"] = filt_sub_tbl["pid"].map(ana_utils.SUB_ID_AGE_DICT)
filt_sub_tbl["sex"] = filt_sub_tbl["pid"].map(ana_utils.SUB_ID_GEN_DICT)
filt_sub_tbl = pd.merge(
    filt_sub_tbl, on="pid", right=ana_utils.SUB_INFO_TBL[["pid", "logwmh"]]
)
filt_sub_tbl.to_csv(f"{out_dir}/sub_list_beforePSM.csv")
