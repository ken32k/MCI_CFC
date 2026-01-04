import os
import glob
import shutil
import numpy as np
import nibabel as nib
import pandas as pd
from joblib import Parallel, delayed
from nilearn import surface
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.datasets import load_conte69


def plot_SNR(sub_id, func_path):
    tSNR_data = surface.load_surf_data(
        f"{func_path}/volumetric/{sub_id}_space-func_desc-se_tSNR.shape.gii"
    )
    print(len(tSNR_data), flush=True)

    inflated_lh = read_surface(f"{FS_OUT_DIR}/{sub_id}/surf/lh.inflated", itype="fs")
    inflated_rh = read_surface(f"{FS_OUT_DIR}/{sub_id}/surf/rh.inflated", itype="fs")

    plot_hemispheres(
        inflated_lh,
        inflated_rh,
        array_name=tSNR_data,
        size=(800, 200),
        color_bar="bottom",
        zoom=1.25,
        embed_nb=True,
        interactive=False,
        share="both",
        cmap="inferno",
        transparent_bg=False,
        screenshot=True,
        filename=f"{DATA_OUT_DIR}/qc/SNR/{sub_id}.png",
    )


def par_proc_fc(sub_id, fd_mean_thr, fd_max_thr, fd_len_thr):
    func_path = f"{MICA_OUT_DIR}/{sub_id}/func/desc-se_dir-AP_rest_bold"
    fc_file = glob.glob(f"{func_path}/surf/{sub_id}*atlas-{ATLAS}_desc-FC.shape.gii")
    fd_file = f"{func_path}/volumetric/{sub_id}_space-func_desc-se_metric_FD.1D"
    fd_png = (
        f"{func_path}/volumetric/{sub_id}_space-func_desc-se_framewiseDisplacement.png"
    )

    if len(fc_file):
        fc_mat = nib.load(fc_file[0]).darrays[0].data
        fc_mat = fc_mat[-200:, -200:]
        np.fill_diagonal(fc_mat, 1)
        row_sums = np.sum(fc_mat, axis=1)

        fd_data = np.loadtxt(fd_file)
        if len(fd_data):
            fd_mean, fd_max, fd_len = np.mean(fd_data), np.max(fd_data), len(fd_data)
        else:
            fd_mean, fd_max, fd_len = np.nan, np.nan, 0

        print(fd_mean, flush=True)

        if fd_mean < fd_mean_thr and fd_max < fd_max_thr and fd_len > fd_len_thr:
            fc_mat = fc_mat + fc_mat.T
            np.fill_diagonal(fc_mat, 1)

            np.save(
                f"{DATA_OUT_DIR}/fcs_raw/{sub_id}_fc-wb_{ATLAS}.npy",
                fc_mat,
            )

            shutil.copyfile(fd_png, f"/public/home/baishw/WMH_MCI/qc/fmri/{sub_id}.png")

            # tSNR_data = surface.load_surf_data(
            #     f"{func_path}/volumetric/{sub_id}_space-func_desc-se_tSNR.shape.gii"
            # )
            # plot_SNR(sub_id, func_path)

            print(
                "[Info].....  ",
                sub_id,
                "Done",
                # np.min(tSNR_data),
                # np.max(tSNR_data),
                flush=True,
            )
            proc_status = 1
        else:
            print("[Error].... " + sub_id + "zero rows exist", flush=True)
            proc_status = 0
    else:
        print("[Error].... " + sub_id + ":no dir", flush=True)
        fd_mean, fd_max, fd_len = np.nan, np.nan, 0
        proc_status = 0

    return {
        "pid": sub_id,
        "proc_status": proc_status,
        "meanFD": fd_mean,
        "maxFD": fd_max,
        "lenFD": fd_len,
    }


def main():
    print(f"Run {__file__}", flush=True)

    global ATLAS, DATA_OUT_DIR, MICA_OUT_DIR, FS_OUT_DIR
    ATLAS = "schaefer-200"  # or "schaefer-400", modify accordingly

    # Set CSVD paths
    HOME = "/public/home/baishw/WMH_MCI"
    DATA_OUT_DIR = f"{HOME}/data"
    MICA_OUT_DIR = f"{HOME}/csvd_mica_out/micapipe_v0.2.0"
    FS_OUT_DIR = f"{HOME}/csvd_mica_out/freesurfer"
    
    fd_mean_thr, fd_max_thr, fd_len_thr = 0.5, 10, 200
    
    proc_mat_list = sorted(os.listdir(MICA_OUT_DIR))
    print("Total %d subjects to process." % len(proc_mat_list))

    os.system(f"find {DATA_OUT_DIR}/fcs_raw/* -type f -delete")
    results = np.array(
        Parallel(n_jobs=-1)(
            delayed(par_proc_fc)(sub_id, fd_mean_thr, fd_max_thr, fd_len_thr)
            for sub_id in proc_mat_list
        )
    )
    print(results)
    pd.DataFrame(pd.json_normalize(results)).to_csv(
        f"{HOME}/data/proc_fc_log.csv", index=False
    )


if __name__ == "__main__":
    main()
