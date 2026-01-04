import ana_utils
import numpy as np
import os
from nilearn import plotting
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed



def plot_subject_cm(sub_id):
    fig_out_dir=f"{ana_utils.PROJ_HOME}/data/cm-fig"
    # print("plot CM: ", sub_id, flush=True)
    # if not os.path.exists(fig_out_dir):
    #     os.makedirs(fig_out_dir)
    # os.system(f"find {fig_out_dir}/* -type f -delete")

    for lm_lab in ana_utils.LM_LAB:
        for cm in ana_utils.CM_ARR:
            cm_lab, cm_idx, _ = cm
            print(cm_idx, cm_lab)
            try:
                fig, axes = plt.subplots(figsize=(8, 8))
                axes.axis("off")
                cm_mat = np.load(
                        f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm_lab}_{cm_lab}_{ana_utils.ATLAS}.npy"
                    )
                cm_mat = cm_mat[cm_idx] if len(np.shape(cm_mat))==3 else cm_mat
                plotting.plot_matrix(
                    mat=cm_mat,
                    colorbar=False,
                    cmap="RdBu_r",
                )
                plt.title("")
                plt.savefig(
                    f"{fig_out_dir}/{sub_id}-{lm_lab}_{cm_lab}.png",
                    bbox_inches="tight",
                )
            except Exception as e:
                print(cm_lab, e, flush=True)



def plot_subject_fc(sub_id):
    print("plot FC: ", sub_id, flush=True)

    try:
        fig, axes = plt.subplots(figsize=(8, 8))
        plotting.plot_matrix(
            mat=np.load(
                f"{ana_utils.PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
            ),
            colorbar=False,
            cmap="CMRmap",
            vmin=-1,
            vmax=1,
        )
        plt.title("")
        axes.axis("off")
        plt.savefig(
            f"{ana_utils.PROJ_HOME}/data/fc-fig/{sub_id}_fc-wb_{ana_utils.ATLAS}.png",
            bbox_inches="tight",
        )
    except Exception as e:
        print(e)


print(f"Run {__file__}", flush=True)
# Get subject list
sub_cm_list = pd.read_csv(f"{ana_utils.PROJ_HOME}/data/sub_list.csv")["pid"].to_list()
# Run parallel
Parallel(n_jobs=-1)(delayed(plot_subject_cm)(sub_id) for sub_id in sub_cm_list)
Parallel(n_jobs=-1)(delayed(plot_subject_fc)(sub_id) for sub_id in sub_cm_list)
