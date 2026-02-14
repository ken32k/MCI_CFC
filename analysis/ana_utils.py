"""Analysis utilities for the MCI_CFC project."""

from __future__ import annotations

import datetime
import glob
import logging
import os
import pickle
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from networkx import *  # noqa: F401,F403 - legacy re-export relied upon downstream
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from nilearn import plotting
import pingouin as pg
import scipy as sp
from scipy.stats import f_oneway
import statsmodels.api as sm
from tqdm import trange
from joblib import Parallel, delayed
from PIL import Image, ImageFont, ImageDraw
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
from neuroCombat import neuroCombat

# brainspace modules
from brainspace.datasets import load_conte69, load_group_fc, load_parcellation
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels

try:  # Optional dependency used in dominance analysis helpers
    from netneurotools.stats import get_dominance_stats
except ImportError:  # pragma: no cover - handled at runtime if unavailable
    get_dominance_stats = None


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# from netneurotools.stats import get_dominance_stats


global PROJ_HOME, DISGROUPS, SUB_ID_LIST, SUB_ID_DIAG_DICT, WMH_LAB, WMH_LAB4, ATLAS, NNODE, TRIU_IDX, COV_TBL, SUB_INFO_DF


def _resolve_project_home() -> Path:
    env_override = os.environ.get("MCI_CFC_HOME")
    base = Path(env_override or "D:/Shanghaitec/PROJECT/CSVD_neo")
    try:
        return base.expanduser().resolve()
    except OSError:
        logger.warning("Falling back to non-resolved project path: %s", base)
        return base


PROJ_PATH = _resolve_project_home()
PROJ_HOME = str(PROJ_PATH)
DATA_DIR = PROJ_PATH / "data"
PARC_DIR = PROJ_PATH / "parc"

# Node and matrix definition
NNODE = 200
ATLAS = f"schaefer-{NNODE}"
TRIU_IDX = np.triu_indices(NNODE, 1)

N_JOBS_DEFAULT = 6
N_JOBS = int(os.environ.get("MCI_CFC_N_JOBS", N_JOBS_DEFAULT))

sub_data_csv = DATA_DIR / "huashan_renji_merge.csv"
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
    ("fg-wei", 2, "#9A7BB7"),
    ("pt-wei", 2, "#7B8CB7"),
    ("si-wei", 2, "#7B8CB7"),
    ("pl-wei", 2, "#999999"),
]

CM_LABS = [cm[0] for cm in CM_ARR[1:]]
CM_ARR_SHORT = [cm[0][:2].upper() for cm in CM_ARR[1:]]
CM_ARR_COLOR = [cm[2] for cm in CM_ARR[1:]]

# Lesion map (LM)
TRACT_LABS = ["wb", "int", "wmh"]
PAL4 = ["#FA7F6f", "#FFBE7A", "#82B0D2"]

def _read_csv_with_encodings(path: Path, encodings: Iterable[str] = ("gbk", "utf-8", "utf-8-sig"), **kwargs) -> pd.DataFrame:
    last_error: Exception | None = None
    for encoding in encodings:
        try:
            return pd.read_csv(path, encoding=encoding, **kwargs)
        except UnicodeDecodeError as exc:
            last_error = exc
    if last_error:
        raise last_error
    raise FileNotFoundError(path)


def _load_subject_dataframe() -> pd.DataFrame:
    df = _read_csv_with_encodings(sub_data_csv)
    df = df[df["Group"].isin(DISGROUPS)].copy()
    logger.info("Subjects from Excel: %s", len(df))
    return df


def _load_subject_ids(candidate_df: pd.DataFrame) -> List[str]:
    sub_list_file = DATA_DIR / "sub_list.csv"
    if sub_list_file.exists():
        ids = _read_csv_with_encodings(sub_list_file)["pid"].to_list()
        if ids:
            return ids
    logger.warning("Subject list file missing or empty, falling back to dataframe IDs")
    return sorted(set(candidate_df["pid"].to_list()))


def _finalize_subject_info(raw_df: pd.DataFrame, subject_ids: List[str]) -> pd.DataFrame:
    df = raw_df[raw_df["pid"].isin(subject_ids)].copy()
    if "gender" in df.columns:
        df["genderv"] = df["gender"].replace({"M": 0, "F": 1})
    output_csv = DATA_DIR / "filt_sub_info.csv"
    try:
        df.to_csv(output_csv, encoding="gbk", index=False)
    except Exception:
        logger.exception("Failed to write filtered subject info to %s", output_csv)
    return df


_RAW_SUB_INFO = _load_subject_dataframe()
SUB_ID_LIST = _load_subject_ids(_RAW_SUB_INFO)
SUB_INFO_DF = _finalize_subject_info(_RAW_SUB_INFO, SUB_ID_LIST)

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


def add_statistic_annotation(p_vals: Sequence[float]) -> List[str]:
    """Return formatted significance annotations for each p-value."""

    def annotate(p: float) -> str:
        if p < 0.001:
            return "< 0.001***"
        if p < 0.01:
            return f"{p:.3f}**"
        if p < 0.05:
            return f"{p:.3f}*"
        return f"{p:.3f}"

    return [annotate(float(p)) for p in p_vals]


def get_stats_sig(p: float) -> str:
    """Map a p-value onto the conventional star-based significance label."""

    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def sub_mat_correlatrion(fc: np.ndarray, sc: np.ndarray) -> Tuple[float, float]:
    """Compute Pearson correlation between FC and SC/CM matrices.

    Args:
        fc: Functional connectivity values (2D matrix or upper-triangle vector).
        sc: Structural/communication matrix with same geometry as fc.

    Returns:
        Tuple of (r, p) Pearson correlation results.
    """

    fc_arr = np.asarray(fc)
    sc_arr = np.asarray(sc)

    if sc_arr.shape == (NNODE, NNODE):
        idx = np.triu_indices(sc_arr.shape[0], k=1)
        fc_arr = fc_arr[idx]
        sc_arr = sc_arr[idx]

    sc_arr = np.where(np.isinf(sc_arr), np.nan, sc_arr)
    valid_mask = ~np.isnan(sc_arr)
    fc_arr = fc_arr[valid_mask]
    sc_arr = sc_arr[valid_mask]

    if fc_arr.size == 0 or sc_arr.size == 0:
        return np.nan, np.nan

    r, p = sp.stats.pearsonr(fc_arr, sc_arr)
    return r, p


def sub_cfc_lr(fc, cms):
    """
    Calculate the adjusted Rsq and total dominance of functional connectome
    and communication matrices correlations for upper triangle.

    Args:
        fc (numpy.ndarray): 2D array representing functional connectome
        cms (numpy.ndarray): 2D array representing communication matrices

    Returns:
        tuple: total dominance, adjusted Rsq list, number of NaN indices
    """
    if get_dominance_stats is None:
        raise ImportError(
            "netneurotools.stats.get_dominance_stats is required for sub_cfc_lr"
        )

    # If 2D array, then select upper triangle and flatten
    cms[np.isinf(cms) | (cms > np.finfo(np.float64).max)] = 0
    nan_indices = np.isnan(cms).any(axis=0)

    X = cms[:, ~nan_indices].T
    y = fc[~nan_indices]

    # Dominance analysis
    model_metrics, model_r_sq = get_dominance_stats(X, y, use_adjusted_r_sq=False)
    return (
        model_metrics["total_dominance"],
        model_r_sq[(0, 1, 2, 3, 4)],
    )


def pg_anova(anova_df, col_name="rsq", factor="Group"):
    """
    ANOVA
    Pingouin one-way ANOVA without covariates. Only for all-fiber tractography comparison.

    Args:
        anova_df (Pandas df): The dataframe for one-way ANOVA
        col_name (str, optional): The column name to be compared. Defaults to "rsq".
        factor (str, optional): The column name representing the factor. Defaults to "Group".

    Returns:
        list: ANOVA and post-hoc results
        f_statistic
        p_value
        posthoc["pval"]
        posthoc["cohen"]
        posthoc["diff"]
    """
    anova_res = pg.anova(data=anova_df, dv=col_name, between=factor, effsize="np2")
    f_statistic = anova_res["F"].iloc[0]
    p_value = anova_res["p-unc"].iloc[0]
    posthoc = pg.pairwise_gameshowell(
        data=anova_df, dv=col_name, between=factor, effsize="cohen"
    )

    return (
        anova_res,
        posthoc,
        f_statistic,
        p_value,
        posthoc["pval"].to_list(),
        posthoc["cohen"].to_list(),
        posthoc["diff"].to_list(),
    )


def pg_mix_anova(anova_df, col_name="rsq", within="Tractography", between="Group", subject="pid"):
    """
    Mixed ANOVA
    Pingouin mixed ANOVA for INT and WMH tractography comparison.
    Args:
        anova_df (Pandas df): The dataframe for mixed ANOVA
        col_name (str, optional): The column name to be compared. Defaults to "rsq".
        within (str, optional): The column name representing the within-subject factor. Defaults to "Tractography".
        between (str, optional): The column name representing the between-subject factor. Defaults to "Group".
        subject (str, optional): The column name representing the subject. Defaults to "pid".
    Returns:
        list: ANOVA and post-hoc results

    """

    anova_res = pg.mixed_anova(
        data=anova_df,
        dv=col_name,
        within=within,
        between=between,
        subject=subject,
        effsize="ng2",
    ).round(3)
    posthoc_df = pd.DataFrame([])
    for tract_lab in TRACT_LABS[1:]:
        lm_anova_df = anova_df[anova_df["Tractography"] == tract_lab]
        posthoc = pg.pairwise_gameshowell(
            data=lm_anova_df, dv=col_name, between="Group", effsize="cohen"
        ).round(3)
        posthoc["Tractography"] = tract_lab.upper()
        posthoc_df = pd.concat([posthoc_df, posthoc], ignore_index=True)

    return (anova_res, posthoc_df)


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


def get_mean_matrix(cm_lab, tract_lab="wb"):
    """get the mean matrix for sc, fc and lesion map"""
    mean_mat = np.zeros((NNODE, NNODE, len(DISGROUPS)))
    mean_n = np.zeros(len(DISGROUPS))
    for sub_id in SUB_ID_LIST:
        if cm_lab == "fc":
            sub_mat = np.load(f"{PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ATLAS}.npy")
        else:
            sub_mat = np.load(
                f"{PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{tract_lab}_{cm_lab}_{ATLAS}.npy"
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
NET_LABS = {
    "Vis": "VN",
    "SomMot": "SMN",
    "DorsAttn": "DAN",
    "SalVentAttn": "VAN",
    "Limbic": "LN",
    "Cont": "FPC",
    "Default": "DMN",
}


def fdr_correction(pvals, alpha=0.05):
    """
    Perform Benjamini-Hochberg FDR correction on a list of p-values.

    Parameters:
    pvals (list or np.array): List or array of p-values to correct.
    alpha (float): Significance level for FDR correction.

    Returns:
    np.array: Array of corrected p-values.
    """
    _, corrected_pvals, _, _ = multipletests(pvals, alpha=alpha, method="fdr_bh")
    return corrected_pvals
