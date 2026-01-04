# use anaconda python 3.8 base on local machine
# use miniconda python 3.11.5 on cluster

import numpy as np
import pandas as pd
import seaborn as sns
import os, sys, glob, datetime
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import scipy as sp

from matplotlib import colors as mcolors
from nilearn import plotting
import pingouin as pg
import pickle
from scipy.stats import f_oneway
import statsmodels.api as sm
from tqdm import trange
from joblib import Parallel, delayed
from PIL import Image, ImageFont, ImageDraw
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols
from sklearn.linear_model import LinearRegression
from neuroCombat import neuroCombat

# brainspace modules
from brainspace.utils.parcellation import map_to_labels
from brainspace.datasets import load_group_fc, load_parcellation
from brainspace.plotting import plot_hemispheres
from brainspace.datasets import load_conte69
# from netneurotools.stats import get_dominance_stats


global PROJ_HOME, DISGROUPS, SUB_ID_LIST, SUB_ID_DIAG_DICT, WMH_LAB, WMH_LAB4, ATLAS, NNODE, TRIU_IDX, COV_TBL, SUB_INFO_DF

# Project directory
PROJ_HOME = os.path("D:/Shanghaitec/PROJECT/CSVD_neo")

# Node and matrix definition
NNODE = 200
ATLAS = f"schaefer-{NNODE}"
TRIU_IDX = np.triu_indices(NNODE, 1)

N_JOBS = 6

sub_data_csv = os.path.join(PROJ_HOME, "data", "huashan_renji_merge.csv")
diag_col = "Group"

DISGROUPS = ["aMCI", "naMCI", "NCI"]
DISGROUPS_SHORT = [disgroup for disgroup in DISGROUPS]

mean_diff_lab = [
    f"{DISGROUPS_SHORT[0]}-{DISGROUPS_SHORT[1]}",
    f"{DISGROUPS_SHORT[0]}-{DISGROUPS_SHORT[2]}",
    f"{DISGROUPS_SHORT[1]}-{DISGROUPS_SHORT[2]}",
]

CM_ARR = [
    ("sc-wei", 0, "#999999"),
    ("co-wei", 0, "#9A7BB7"),
    ("fg-wei", 3, "#9A7BB7"),
    ("pt-wei", 3, "#7B8CB7"),
    ("si-wei", 3, "#7B8CB7"),
    ("pl-wei", 3, "#999999"),
]

CM_LAB_ARR = [cm[0] for cm in CM_ARR[1:]]
CM_ARR_SHORT = [cm[0][:2].upper() for cm in CM_ARR[1:]]
CM_ARR_COLOR = [cm[2] for cm in CM_ARR[1:]]

# Lesion map (LM)
LM_LAB = ["wb", "int", "wmh"]
PAL4 = ["#FA7F6f", "#FFBE7A", "#82B0D2"]

# Clinical information
SUB_INFO_DF = pd.read_csv(sub_data_csv, encoding="gbk")
SUB_INFO_DF = SUB_INFO_DF[SUB_INFO_DF["Group"].isin(DISGROUPS)]
print("Subjects from Excel: ", len(SUB_INFO_DF))

# Get subject list
sub_list_file = f"{PROJ_HOME}/data/sub_list.csv"
if os.path.exists(sub_list_file):
    SUB_ID_LIST = pd.read_csv(sub_list_file, encoding="gbk")["pid"].to_list()
else:
    SUB_ID_LIST = list(set(SUB_INFO_DF.pid.tolist()))

# Filter SUB_INFO_DF
SUB_INFO_DF = SUB_INFO_DF[SUB_INFO_DF["pid"].isin(SUB_ID_LIST)]
SUB_INFO_DF["genderv"] = SUB_INFO_DF["gender"].replace({"M": 0, "F": 1})
SUB_INFO_DF.to_csv(f"{PROJ_HOME}/data/filt_sub_info.csv", encoding="gbk")

# Get the diagnosis dictionary
SUB_ID_DIAG_DICT = SUB_INFO_DF.set_index("pid")["Group"].to_dict()
SUB_ID_AGE_DICT = SUB_INFO_DF.set_index("pid")["age"].to_dict()
SUB_ID_GEN_DICT = SUB_INFO_DF.set_index("pid")["genderv"].to_dict()
SUB_ID_WMH_DICT = SUB_INFO_DF.set_index("pid")["logwmh"].to_dict()
COV_TBL = SUB_INFO_DF[["pid", "age", "genderv", "logwmh"]]

value_counts = {
    value: sum(1 for v in SUB_ID_DIAG_DICT.values() if v == value)
    for value in set(SUB_ID_DIAG_DICT.values())
}
print(
    "SUB_ID_LIST: list, format is like this: ",
    SUB_ID_LIST[0],
    ", and total length of SUB_ID_LIST =",
    len(SUB_ID_LIST),
    value_counts,
)


def add_statistic_annotation(p_vals):
    """
        Add statistical annotation based on p-values.
        Args:
            p_vals (float or list): A single p-value or a list of p-values.
        Returns:
            list or str: Statistical annotations as a list (if input is a list) or a string (if input is a float).
    """
    def annotate(p):
        if p < 0.001:
            return "< 0.001***"
        elif p < 0.01:
            return f"{p}**"
        elif p < 0.05:
            return f"{p}*"
        else:
            return f"{p}"

    if isinstance(p_vals, list):
        return [annotate(p) for p in p_vals]
    elif isinstance(p_vals, (float, int)):
        return annotate(p_vals)
    else:
        raise ValueError("p_vals must be a float, int, or a list of floats/ints.")


def sub_mat_correlatrion(fc, sc):
    """
    Calculate the rsq of fc and matrices correlations
    FC should be inputed before SC!!!
    """

    # If 2D array, then select upper triangle and flatten
    if sc.shape == (NNODE, NNODE):
        fc, sc = (
            fc[np.triu_indices(fc.shape[0], k=1)],
            sc[np.triu_indices(sc.shape[0], k=1)],
        )

    # Replace inf and 0 values with nan in sc
    sc = np.where(np.isinf(sc), 0, sc)

    # Remove Nans according to nans in sc
    fc, sc = fc[~np.isnan(sc)], sc[~np.isnan(sc)]

    # Z score (unnecessary)
    # fc = sp.stats.zscore(fc)
    # sc = sp.stats.zscore(sc)

    # pearson correlation
    r, p = sp.stats.pearsonr(fc, sc)
    # spearman correlation (not recommended)
    # r, p = sp.stats.spearmanr(fc, sc)

    return r, p


def sub_mat_linear_reg(fc, cms):
    """
    Calculate the adjusted Rsq and total dominance of functional connectome
    and communication matrices correlations for upper triangle.

    Args:
        fc (numpy.ndarray): 2D array representing functional connectome
        cms (numpy.ndarray): 2D array representing communication matrices

    Returns:
        list: Adjusted Rsq, total dominance
    """
    # Linear regression
    cms[np.isinf(cms) | (cms > np.finfo(np.float64).max)] = 0
    nan_indices = np.isnan(cms).any(axis=0)

    X = cms[:, ~nan_indices].T
    y = fc[~nan_indices]

    # Dominance analysis
    model_metrics, model_r_sq = get_dominance_stats(X, y, use_adjusted_r_sq=False)
    return (
        0,
        0,
        model_metrics["total_dominance"],
        model_r_sq[(0, 1, 2, 3, 4)],
        len(nan_indices),
    )


def pg_anova(anova_df, col_name="rsq", factors="Group"):
    """
    ANOVA
    Pingouin one-way ANOVA without covariates. Only for all-fiber tractography comparison.

    Args:
        anova_df (Pandas df): The dataframe for one-way ANOVA
        col_name (str, optional): The column name to be compared. Defaults to "rsq".
        factors (str, optional): The column name representing the factor. Defaults to "Group".

    Returns:
        list: ANOVA and post-hoc results
        f_statistic
        p_value
        posthoc["pval"]
        posthoc["cohen"]
        posthoc["diff"]
    """
    anova_df[factors] = pd.Categorical(anova_df[factors], categories=DISGROUPS, ordered=True)
    anova_res = pg.anova(data=anova_df, dv=col_name, between=factors, effsize="np2").round(3)
    f_statistic = anova_res["F"].iloc[0]
    p_value = anova_res["p-unc"].iloc[0]
    posthoc = pg.pairwise_gameshowell(data=anova_df, dv=col_name, between=factors, effsize="cohen").round(3)

    return (
        anova_res,
        posthoc,
        f_statistic,
        p_value,
        posthoc["pval"].to_list(),
        posthoc["cohen"].to_list(),
        posthoc["diff"].to_list(),
    )
    

def pg_mix_anova(anova_df, col_name="rsq", within="lm", between="Group", subject="pid"):
    """
    Mixed ANOVA
    Pingouin mixed ANOVA for INT and WMH tractography comparison.
    Args:
        anova_df (Pandas df): The dataframe for mixed ANOVA
        col_name (str, optional): The column name to be compared. Defaults to "rsq".
        within (str, optional): The column name representing the within-subject factor. Defaults to "lm".
        between (str, optional): The column name representing the between-subject factor. Defaults to "Group".
        subject (str, optional): The column name representing the subject. Defaults to "pid".
    Returns:
        list: ANOVA and post-hoc results
        f_statistic
        p_value
        posthoc["pval"]

    """

    anova_res = pg.mixed_anova(data=anova_df, dv=col_name, within=within, between=between,
                                subject=subject, effsize="ng2").round(3)
    f_statistic = anova_res["F"].iloc
    p_value = anova_res["p-unc"].iloc
    # Check for significant interaction of Group 
    posthoc = None
    if p_value[0] < 0.05:
        print("Significant Group effect detected.")
        
        for lm_lab in LM_LAB[1:]:
            lm_anova_df = anova_df[anova_df['lm'] == lm_lab]
            posthoc = pg.pairwise_gameshowell(data=lm_anova_df, dv='rsq', between='Group', effsize='cohen').round(3)
            print(f"\nPost-hoc comparisons (averaged across {lm_lab.upper()}-tractography):")
            display(posthoc)
            
    # Check for significant main effect of Tractography
    # if p_value[1] < 0.05:
    #     print("Significant main effect of Tractography detected.")
    #     for group in anova_df['Group'].unique():
    #         df_group = anova_df[anova_df['Group'] == group]
    #         pivot_df = df_group.pivot(index='pid', columns='lm', values=col_name)
            
    #         # Paired t-test between INT and WMH within each group
    #         ttest_result = pg.ttest(pivot_df['int'], pivot_df['wmh'], paired=True, correction='auto',
    #                     alternative='two-sided').round(3)
    #         display(ttest_result)
            
    if p_value[2] < 0.05:
        print("Significant interaction between LM and Group detected detected.")
        posthoc = pg.pairwise_gameshowell(data=anova_df, dv=col_name, between=between, effsize="cohen").round(3)
        print("\nPost-hoc comparisons:")
        display(posthoc)
    return (
        anova_res,
        posthoc,
        f_statistic,
        p_value,
        posthoc["pval"].to_list() if posthoc is not None else [],
        posthoc["cohen"].to_list() if posthoc is not None else [],
        posthoc["diff"].to_list() if posthoc is not None else [],
    )
    
    
def get_yeo7_dict(order_file):
    
    """
    Get the yeo 7 or 17 subnetwork-node map

    Args:
        order_file (file): _description_

    Returns:
        yeo7_dict (diction): _description_
    """
    yeo7_dict = dict()
    order_tbl = pd.read_csv(order_file, header=None, sep="\t")[[0, 1]]
    order_tbl.columns = ["idx", "node"]
    order_tbl["net"] = order_tbl["node"].apply(lambda x: x.split("_")[2])
    order_tbl["idx"] = order_tbl["idx"] - 1
    yeo7_dict = order_tbl.groupby("net")["idx"].apply(list).to_dict()
    return yeo7_dict

def get_mean_matrix(cm_lab, lm="wb"):
    """get the mean matrix for sc, fc and lesion map"""
    mean_mat = np.zeros((NNODE, NNODE, len(DISGROUPS)))
    mean_n = np.zeros(len(DISGROUPS))
    for sub_id in SUB_ID_LIST:
        if cm_lab == "fc":
            sub_mat = np.load(f"{PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ATLAS}.npy")
        else:
            sub_mat = np.load(
                f"{PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm}_{cm_lab}_{ATLAS}.npy"
            )[0]

        mean_idx = DISGROUPS.index(SUB_ID_DIAG_DICT[sub_id])
        mean_mat[:, :, mean_idx] += sub_mat
        mean_n[mean_idx] += 1
    mean_n = np.where(mean_n != 0, mean_n, 1)
    mean_mat = mean_mat / mean_n[None, None, :]
    return (mean_mat, mean_n)


def plot_single_matrix(mat, output_path, title, cmap="BuPu", vmin=0, vmax=50):
    fig, ax = plt.subplots()
    im = ax.imshow(mat, cmap=cmap, vmax=vmax, vmin=vmin)
    plt.axis("off")
    plt.title(title)
    plt.savefig(output_path, bbox_inches="tight", dpi=300)


def plot_ax_swarmplot(
    plot_data, ax, i=0, x="Group", y="rsq", ylim=(0, 0.2), boxwidth=0.4, swarmsize=1.5
):
    """
    plot_ax_swarmplot and box

    Args:
        plot_data (obj): pandas DataFrame
        ax (_type_): plt ax
        i (int, optional): _description_. Defaults to 0.
        x (str, optional): _description_. Defaults to "Group".
        y (str, optional): _description_. Defaults to "rsq".
        ylim (turple, optional): _description_. Defaults to (0, 0.2).
        swarmsize (int, optional): _description_. Defaults to 2.

    Returns:
        ax: plt ax
    """
    sns.swarmplot(
        data=plot_data,
        x="Group",
        y=y,
        ax=ax,
        size=swarmsize,
        order=DISGROUPS,
        color="gray",
        alpha=0.9,
        palette=PAL4,
        zorder=1,
    )
    sns.boxplot(
        data=plot_data,
        x="Group",
        y=y,
        linecolor="k",
        ax=ax,
        width=boxwidth,
        order=DISGROUPS,
        fliersize=0,
        palette=PAL4,
        zorder=2,
    )

    if i > 12:
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.set_ylabel("")
    ax.set_ylim(ylim)
    ax.set_xlabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels(DISGROUPS_SHORT, rotation=90)
    ax.set_facecolor("whitesmoke")

    return ax


def plot_create_image_grid(image_files, nrow, ncol, padding, output_path, is_dlt=False):
    """
    Concat multiple image files

    Args:
        image_files (string): image files path
        nrow (int): number of rows
        ncol (int): number of columns
        padding (int): padding in pixels
        output_path (string): output path
        is_dlt (bool, optional): Delete after merge. Defaults to False.
    """

    img_list = [Image.open(f) for f in image_files]
    width, height = img_list[0].size
    final_width = ncol * width + (ncol) * 2 * padding
    final_height = nrow * height + (nrow) * 2 * padding
    new_im = Image.new("RGB", (final_width, final_height), color="white")
    for index, im in enumerate(img_list):
        x_offset = (index % ncol) * (width + padding)
        y_offset = (index // ncol) * (height + padding)
        new_im.paste(im, (x_offset, y_offset))
    new_im.save(output_path)
    if is_dlt:
        for file in image_files:
            os.remove(file)


YEO7_DICT = get_yeo7_dict(
    f"{PROJ_HOME}/parc/Schaefer2018_200Parcels_7Networks_order.txt",
)
net_abbr = {
    "Vis": "VN",
    "SomMot": "SMN",
    "DorsAttn": "DAN",
    "SalVentAttn": "VAN",
    "Limbic": "LN",
    "Cont": "FPC",
    "Default": "DMN",
}

