import ana_utils
import numpy as np
import pandas as pd

import os, sys
import scipy as sp
import nilearn as nil
from nilearn import plotting
from scipy.stats import f_oneway
import statsmodels.api as sm

import matplotlib.pyplot as plt
import seaborn as sns


# Filter subjects by sub_list
sub_id_list = ana_utils.SUB_ID_LIST
SUB_INFO_TBL = ana_utils.SUB_INFO_TBL
mean_diff_lab = ana_utils.mean_diff_lab

output_dir = f"{ana_utils.PROJ_HOME}/results/lms_stat"

nawm_tbl = pd.read_csv(
    f"{ana_utils.PROJ_HOME}/results/lms_stat/nawm_volume_raw.csv", names=["pid", "nawm"]
)
nawm_tbl["lognawm"] = nawm_tbl["nawm"].apply(lambda x: np.log(x))
SUB_INFO_TBL = pd.merge(SUB_INFO_TBL, on="pid", right=nawm_tbl)

wmh_metric_tbl = pd.read_csv(
    f"{ana_utils.PROJ_HOME}/results/lms_stat/4_wmh_dwi_metric.csv",
    names=["pid", "dwim", "val"],
)
wmh_metric_tbl = wmh_metric_tbl[wmh_metric_tbl.val > 0]
wmh_metric_pivot = wmh_metric_tbl.pivot_table(
    index="pid", values="val", columns="dwim"
).reset_index()
SUB_INFO_TBL = pd.merge(SUB_INFO_TBL, on="pid", right=wmh_metric_pivot)


SUB_INFO_TBL["Group"] = SUB_INFO_TBL["pid"].map(ana_utils.SUB_ID_DIAG_DICT)
SUB_INFO_TBL = SUB_INFO_TBL[SUB_INFO_TBL.pid.isin(sub_id_list)]
SUB_INFO_TBL = SUB_INFO_TBL.dropna(subset=["Group"])
comp_res = []



# Plot
plot_items = ["age", "logwmh", "lognawm", "fa", "rd", "adc",]
plot_titles = [
    "Age",
    "log(WMH Volume)",
    "log(NAWM Volume)",
    "FA(WMH)",
    "RD(WMH)",
    "MD(WMH)",
]

plot_ylims = [(40, 100), (1.5, 6.5), (12, 14), (0.10, 0.60), (0.0, 0.003), (0, 0.003)]

fig, axes = plt.subplots(2, 3, figsize=(6,4))
for i, plot_item in enumerate(plot_items):

    anova_f, anova_p, group_p, _, posthoc_p, meandiffs  = ana_utils.pg_anova(SUB_INFO_TBL, plot_item)
    print(posthoc_p)
    comp_res.append(
        {
            "cm": "0",
            "lm": plot_item,
            "anova_f": anova_f,
            "group_p": group_p,
            "posthoc_p_0": posthoc_p[0],
            "posthoc_p_1": posthoc_p[1],
            "posthoc_p_2": posthoc_p[2],
            mean_diff_lab[0]: meandiffs[0],
            mean_diff_lab[1]: meandiffs[1],
            mean_diff_lab[2]: meandiffs[2],
        }
    )

    ax = axes.flatten()[i]
    ax= ana_utils.plot_ax_swarmplot(SUB_INFO_TBL, ax, i=0, x=ana_utils.diag_col, y=plot_item, ylim=plot_ylims[i])
    sub_item_data = SUB_INFO_TBL[plot_item]
    ylim_top = np.max(sub_item_data) + 3 * np.std(sub_item_data)
    ylim_bottom = np.min(sub_item_data) - 3 * np.std(sub_item_data)
    ylim_bottom, ylim_top = plot_ylims[i]
    ax.set_yticks(np.round(np.linspace(ylim_bottom, ylim_top, 5), 4))
    ax.set_xlabel("")
    ax.set_ylabel(plot_titles[i])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels(ana_utils.DISGROUPS_SHORT, rotation=90)
    ax.set_facecolor("whitesmoke")
    # Add significance
    if ax.get_ylim()[0] > ax.get_ylim()[1]:
        ax.invert_yaxis()
    sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
    for i in range(len(posthoc_p)):
        if group_p < 0.05 and posthoc_p[i] < 0.05:
            ax.text(
                0.5 + 0.5 * i,
                ylim_top * (0.95 - i * 0.05),
                ana_utils.get_stats_sig(group_p),
                ha="center",
            )
            ax.axhline(
                ylim_top * (0.95 - i * 0.05),
                xmin=sig_line_xs[i][0],
                xmax=sig_line_xs[i][1],
                color="black",
                linestyle="-",
                zorder=1,
            )
    
plt.tight_layout()
plt.show()
plt.savefig(
    f"{output_dir}/wmh-vol_plot.svg",
    bbox_inches="tight",
)

# ---------------------

# Filter subjects by sub_list
sub_id_list = pd.read_csv(f"{ana_utils.PROJ_HOME}/data/sub_list.csv")["pid"].to_list()

output_dir = f"{ana_utils.PROJ_HOME}/results/lms_stat"
SUB_INFO_TBL = ana_utils.SUB_INFO_TBL

nawm_tbl = pd.read_csv(
    f"{ana_utils.PROJ_HOME}/results/lms_stat/nawm_volume_raw.csv", names=["pid", "nawm"]
)
nawm_tbl["lognawm"] = nawm_tbl["nawm"].apply(lambda x: np.log(x))
SUB_INFO_TBL = pd.merge(SUB_INFO_TBL, on="pid", right=nawm_tbl)

wmh_metric_tbl = pd.read_csv(
    f"{ana_utils.PROJ_HOME}/results/lms_stat/2_wmh_dwi_metric.csv",
    names=["pid", "dwim", "val"],
)
wmh_metric_tbl = wmh_metric_tbl[wmh_metric_tbl.val > 0]
wmh_metric_pivot = wmh_metric_tbl.pivot_table(
    index="pid", values="val", columns="dwim"
).reset_index()
SUB_INFO_TBL = pd.merge(SUB_INFO_TBL, on="pid", right=wmh_metric_pivot)


SUB_INFO_TBL["Group"] = SUB_INFO_TBL["pid"].map(ana_utils.SUB_ID_DIAG_DICT)
SUB_INFO_TBL = SUB_INFO_TBL[SUB_INFO_TBL.pid.isin(sub_id_list)]
SUB_INFO_TBL = SUB_INFO_TBL.dropna(subset=["Group"])

# Plot
plot_items = ["age", "logwmh", "lognawm", "fa", "rd", "adc"]
plot_titles = [
    "Age",
    "log(WMH Volume)",
    "log(NAWM Volume)",
    "FA(NAWM)",
    "RD(NAWM)",
    "MD(NAWM)",
]
plot_ylims = [(40, 100), (1.5, 6.5), (12, 14), (0.10, 0.6), (0.0005, 0.001), (0.0007, 0.0012)]
fig, axes = plt.subplots(2, 3, figsize=(6,4))
for i, plot_item in enumerate(plot_items):
    
    anova_f, anova_p, group_p, _, posthoc_p, meandiffs  = ana_utils.pg_anova(SUB_INFO_TBL, plot_item)
    comp_res.append(
        {
            "cm": "0",
            "lm": plot_item,
            "anova_f": anova_f,
            "group_p": group_p,
            "posthoc_p_0": posthoc_p[0],
            "posthoc_p_1": posthoc_p[1],
            "posthoc_p_2": posthoc_p[2],
            mean_diff_lab[0]: meandiffs[0],
            mean_diff_lab[1]: meandiffs[1],
            mean_diff_lab[2]: meandiffs[2],
        }
    )

    ax = axes.flatten()[i]
    ax= ana_utils.plot_ax_swarmplot(SUB_INFO_TBL, ax, i=0, x=ana_utils.diag_col, y=plot_item, ylim=plot_ylims[i])


    sub_item_data = SUB_INFO_TBL[plot_item]
    ylim_top = np.max(sub_item_data) + 3 * np.std(sub_item_data)
    ylim_bottom = np.min(sub_item_data) - 3 * np.std(sub_item_data)
    ylim_bottom, ylim_top = plot_ylims[i]
    ax.set_yticks(np.round(np.linspace(ylim_bottom, ylim_top, 5), 4))
    ax.set_xlabel("")
    ax.set_ylabel(plot_titles[i])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels(ana_utils.DISGROUPS_SHORT, rotation=90)
    ax.set_facecolor("whitesmoke")
    # Add significance
    if ax.get_ylim()[0] > ax.get_ylim()[1]:
        ax.invert_yaxis()
    sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
    for i in range(len(posthoc_p)):
        if group_p < 0.05 and posthoc_p[i] < 0.05:
            ax.text(
                0.5 + 0.5 * i,
                ylim_top * (0.95 - i * 0.05),
                ana_utils.get_stats_sig(group_p),
                ha="center",
            )
            ax.axhline(
                ylim_top * (0.95 - i * 0.05),
                xmin=sig_line_xs[i][0],
                xmax=sig_line_xs[i][1],
                color="black",
                linestyle="-",
                zorder=1,
            )
    
plt.tight_layout()
plt.show()
plt.savefig(
    f"{output_dir}/nawm-vol_plot.svg",
    bbox_inches="tight",
)


pd.DataFrame(comp_res).to_csv(output_dir+"/_comp_results.csv")