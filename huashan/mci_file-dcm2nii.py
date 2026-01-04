"""
Author: Ken32g ken32k@163.com
Date: 2024-08-13 16:50:38
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2024-08-13 16:52:54
FilePath: /csvd-sfc/huashan/mci_file-dcm2nii.py
Description: This file transform the dcm files to nifti file, only for huashan datasets.
"""

import os, glob
import pandas as pd
import subprocess
import fnmatch
from joblib import Parallel, delayed


PROJ_HOME = "/public_bme/data/v-baishw/CSVD"
# PROJ_HOME = "/public/home/baishw/WMH_MCI"
bids_dir = f"{PROJ_HOME}/csvd_mica_bids"

sub_list_path = f"{PROJ_HOME}/data/sub_list_huashan_raw.csv"
# sub_list_path = f"{PROJ_HOME}/data/sub_list_psm2515.csv"

huashan_data_path = "/public_bme/share/HUASHAN/PET_center/MRI"
# huashan_data_path = "/public_bme/share/HUASHAN/PET_center/MRI/data"

# print(len(SUB_INFO_TBL))
sub_list_psm = pd.read_csv(sub_list_path, encoding="utf-8")
# sub_list_psm = sub_list_psm[sub_list_psm.hospital == "huashan"]
SUB_INFO_TBL = sub_list_psm.copy()
# SUB_INFO_TBL = SUB_INFO_TBL.head(5)
print(len(SUB_INFO_TBL), flush=True)

seq_keyword, scan_dir, des_bid_dir, des_seq_name = "dark", "_", "anat", "FLAIR"

# seq_keyword, scan_dir, des_bid_dir, des_seq_name = "mprage", "_", "anat", "T1w"

# seq_keyword, scan_dir, des_bid_dir, des_seq_name = (
#     "BOLD",
#     "_PA_",
#     "func",
#     "dir-PA_rest_bold",
# )
# seq_keyword, scan_dir, des_bid_dir, des_seq_name = (
#     "BOLD",
#     "_AP_",
#     "func",
#     "dir-AP_rest_bold",
# )
# seq_keyword, scan_dir, des_bid_dir, des_seq_name = (
#     "EP2D_DIFF",
#     "_AP_",
#     "dwi",
#     "dir-AP_dwi",
# )
# seq_keyword, scan_dir, des_bid_dir, des_seq_name = (
#     "EP2D_DIFF",
#     "_PA_",
#     "dwi",
#     "dir-PA_dwi",
# )


def dcm2nii(row):
    sub_pid, sub_path_name = row["pid"], row["path_name"]
    if os.path.exists(
        f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}/{des_bid_dir}/{sub_pid}_{des_seq_name}.nii.gz"
    ):
        print("[Info] ... Skip:", sub_pid, flush=True)
        return

    for sub_dir in os.listdir(huashan_data_path):
        target_dir_exist = True
        # Check if sub_path_name is part of sub_dir
        if sub_path_name in sub_dir:
            target_dir_exist = 1
            print("[Info] ... Found:", sub_pid, sub_path_name, flush=True)
            # Construct the full path to the subdirectory
            full_path = os.path.join(huashan_data_path, sub_dir)
            # Walk through the directory
            for dirpath, _, _ in os.walk(full_path):
                if (
                    f"{seq_keyword}".lower() in dirpath.lower()
                    and f"{scan_dir}".lower() in dirpath.lower()
                    and os.path.getsize(dirpath)
                ):

                    target_dir_exist = True
                    subprocess.run(
                        ["mkdir", "-p", f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}"]
                    )

                    subprocess.run(
                        [
                            "mkdir",
                            "-p",
                            f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}/{des_bid_dir}",
                        ]
                    )
                    subprocess.run(
                        ["rm", "-rf", f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp"]
                    )
                    subprocess.run(
                        ["mkdir", "-p", f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp"]
                    )
                    try:
                        subprocess.run(
                            [
                                "dcm2niix",
                                "-f",
                                "%d",
                                "-z",
                                "o",
                                "-o",
                                f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp/",
                                f"{dirpath}",
                                # ">/dev/null",
                            ]
                        )

                        cmd_mv_nii = f"mv {PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp/*.nii.gz {PROJ_HOME}/csvd_mica_bids/{sub_pid}/{des_bid_dir}/{sub_pid}_{des_seq_name}.nii.gz"
                        subprocess.run(cmd_mv_nii, shell=True, check=True)
                        cmd_mv_json = f"mv {PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp/*.json {PROJ_HOME}/csvd_mica_bids/{sub_pid}/{des_bid_dir}/{sub_pid}_{des_seq_name}.json"
                        subprocess.run(cmd_mv_json, shell=True, check=True)
                        # Move bval and bvec
                        if des_bid_dir == "dwi":
                            cmd_mv_json = f"mv {PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp/*.bval {PROJ_HOME}/csvd_mica_bids/{sub_pid}/{des_bid_dir}/{sub_pid}_{des_seq_name}.bval"
                            subprocess.run(cmd_mv_json, shell=True, check=True)
                            cmd_mv_json = f"mv {PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp/*.bvec {PROJ_HOME}/csvd_mica_bids/{sub_pid}/{des_bid_dir}/{sub_pid}_{des_seq_name}.bvec"
                            subprocess.run(cmd_mv_json, shell=True, check=True)

                        subprocess.run(
                            ["rm", "-rf", f"{PROJ_HOME}/csvd_mica_bids/{sub_pid}/tmp"]
                        )

                    except:
                        print("err")
                    break
        else:
            target_dir_exist = False

    if not target_dir_exist:
        print("[Info] ... Not Found:", sub_pid, sub_path_name, flush=True)


n = 0
Parallel(n_jobs=-1)(delayed(dcm2nii)(row) for _, row in SUB_INFO_TBL.iterrows())


print(os.listdir("/public_bme/data/v-baishw/CSVD/csvd_mica_bids/sub-642/dwi/"))
