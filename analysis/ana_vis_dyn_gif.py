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
import subprocess

from PIL import Image

cm_arr = ["fg-wei"]
for cm_lab in cm_arr:
    for sub_id in [ "sub-376", "sub-614", "sub-137","sub-283"]:
    # 
        fc_mat = np.load(
            f"{ana_utils.PROJ_HOME}/data/fcs/" + sub_id + f"_fc-wb_{ana_utils.ATLAS}.npy"
        )

        # Plot the FC and save
        ana_utils.plot_single_matrix(
            fc_mat,
            output_path=f"{ana_utils.PROJ_HOME}/results/single/{cm_lab}/{sub_id}_fc.png",
            title="Subject: %s, Functional Connectivity" % (str(sub_id)),
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
        )

        cm_mats = np.load(
            f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_wb_{cm_lab}_{ana_utils.ATLAS}.npy"
        )

        # if "fg" in cm_lab:
        #     markov_t = np.arange(0.2, 15.2, 0.2)
        markov_t = np.arange(0.2, 15.2, 0.2)
        images = []

        for i, frame in enumerate(cm_mats):
            frame = (frame * 255).astype(np.uint8)

            ana_utils.plot_single_matrix(
                frame,
                output_path=f"{ana_utils.PROJ_HOME}/results/single/{cm_lab}/tmp-{sub_id}-{str(i)}.png",
                title="Subject: %s, markov time = %.2f" % (sub_id, markov_t[i]),
            )
            images.append(
                Image.open(
                    f"{ana_utils.PROJ_HOME}/results/single/{cm_lab}/tmp-{sub_id}-{str(i)}.png"
                )
            )

        # Delete tmp png files
        # subprocess.run(f"rm {ana_utils.PROJ_HOME}/results/single/{cm_lab}/tmp-*.png", shell=True)

        # Save GIF
        images[0].save(
            f"{ana_utils.PROJ_HOME}/results/single/{cm_lab}/{sub_id}_{cm_lab}.gif",
            save_all=True,
            append_images=images[1:],
            duration=200,
            loop=0,
        )
