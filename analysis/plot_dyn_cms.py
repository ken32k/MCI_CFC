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


markov_ts = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
gammas = [0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
out_dir = f"{ana_utils.PROJ_HOME}/results/corr_glob_lm_cfc"

cm_arr = [
    "fg-wei",
    "pl-wei",
    "si-wei",
    "pt-wei",
]
for lm_lab in ana_utils.LM_LAB:
    for cm_lab in cm_arr:
        print(cm_lab)

        glob_cmfc_tbl = pd.read_csv(f"{out_dir}/{cm_lab}-cfc_{lm_lab}.csv")
        glob_cmfc_tbl["significance"] = glob_cmfc_tbl['pval'].apply(lambda x: 1 if x < 0.05 else 0)
        
        plt.figure(figsize=(12, 6))
        ax = plt.gca()
        sns.countplot(data=glob_cmfc_tbl,
            x="idx",hue="significance", stat="percent",
            palette="Spectral")
        plt.show()
        plt.savefig(
            f"{out_dir}/{cm_lab}-sfc_{lm_lab}_pvals.png",
            bbox_inches="tight",
        )
        
        # Lineplot overall
        plt.figure(figsize=(12, 6))
        ax = plt.gca()
        sns.lineplot(
            data=glob_cmfc_tbl,
            x="idx",
            y="rsq",
            errorbar="sd",
            palette="Spectral",
        )
        # ax.set_xscale("log")

        # add title
        if "fg" in cm_lab:
            plt.title(f"{cm_lab}-FC Coupling by t")
            plt.xlabel("t")
            plt.ylabel("R")
            # plt.xticks(markov_ts, labels=markov_ts)
        else:
            plt.title(f"{cm_lab}-FC Coupling by gamma")
            plt.xlabel("gamma")
            plt.ylabel("R")
            # plt.xticks(gammas, labels=gammas)

        plt.show()
        plt.savefig(
            f"{out_dir}/{cm_lab}-sfc_{lm_lab}.png",
            bbox_inches="tight",
        )

        
        # Lineplot of subgroups
        plt.figure(figsize=(12, 6))
        ax = plt.gca()
        sns.lineplot(
            data=glob_cmfc_tbl,
            x="idx",
            y="rsq",
            hue="Group",
            errorbar="sd",
            hue_order=ana_utils.DISGROUPS,
            palette="Spectral",
        )
        # ax.set_xscale("log")

        # add title
        if "fg" in cm_lab:
            plt.title(f"{cm_lab}-FC Coupling by t")
            plt.xlabel("t")
            plt.ylabel("R")
            # plt.xticks(markov_ts, labels=markov_ts)
        else:
            plt.title(f"{cm_lab}-FC Coupling by gamma")
            plt.xlabel("gamma")
            plt.ylabel("R")
            # plt.xticks(gammas, labels=gammas)

        plt.show()
        plt.savefig(
            f"{out_dir}/{cm_lab}-sfc_{lm_lab}.png",
            bbox_inches="tight",
        )
