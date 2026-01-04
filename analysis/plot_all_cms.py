import ana_utils
import numpy as np
import pandas as pd

import os, sys
import scipy as sp
import nilearn as nil
from nilearn import plotting
from scipy.stats import f_oneway
import statsmodels.api as sm
from matplotlib.ticker import StrMethodFormatter
import matplotlib.pyplot as plt
import seaborn as sns

sfc_method = "corr_glob"
out_dir = f"{ana_utils.PROJ_HOME}/results/{sfc_method}_lm_cfc"
plot_ylim = 0.24

mean_diff_lab = ana_utils.mean_diff_lab
cm_arr = ana_utils.CM_ARR[1:]
lm_labs = ana_utils.LM_LAB
comp_res = []


# Plot all cfcs
for lm_idx, lm_lab in enumerate(lm_labs):
    print(f"[Info] Plot {lm_lab}", flush=True)
    num_subplots = len(cm_arr)
    num_columns = 8

    num_rows = num_subplots // num_columns + 1
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(8, 2.5))
    axes = axes.flatten()
    plt.subplots_adjust(wspace=0.2)
    for i in range(num_subplots):
        cm_lab, cm_idx, _ = cm_arr[i]

        glob_cmfc_tbl = pd.read_csv(f"{out_dir}/{cm_lab}-cfc_{lm_lab}.csv")
        glob_cmfc_tbl = glob_cmfc_tbl[glob_cmfc_tbl["idx"] == cm_idx]
        # ANOVA
        anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_anova(
            glob_cmfc_tbl
        )
        # Append ANOVA results
        comp_res.append(
            {
                "cm": cm_lab,
                "lm": lm_lab,
                "anova_f": anova_f,
                "anova_p": anova_p,
                "p_group": p_group,
                "p_age": p_age,
                "posthoc_p_0": posthoc_p[0],
                "posthoc_p_1": posthoc_p[1],
                "posthoc_p_2": posthoc_p[2],
                mean_diff_lab[0]: meandiffs[0],
                mean_diff_lab[1]: meandiffs[1],
                mean_diff_lab[2]: meandiffs[2],
            }
        )
        ax = axes[i]
        ax = ana_utils.plot_ax_swarmplot(
            glob_cmfc_tbl, ax, i=0, x="Group", y="rsq", ylim=plot_ylim
        )

        ax.set_ylabel(["All", "INT", "WMH"][lm_idx])
        ax.set_ylim(0, plot_ylim)
        ax.set_yticks(np.round(np.linspace(0, plot_ylim, 4), 3))
        sns.despine(trim=True, ax=ax, offset=1)
        if i > 0:
            ax.spines["left"].set_visible(False)
            ax.set(yticks=[], ylabel="")

        ax.set_title(cm_lab[:2].upper())

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("")
        ax.set_xticklabels(ana_utils.DISGROUPS_SHORT, rotation=90)

        sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
        for i in range(len(posthoc_p)):
            if p_group < 0.05 and posthoc_p[i] < 0.05:
                ax.text(
                    0.5 + 0.5 * i,
                    plot_ylim * (0.9 - i * 0.1),
                    ana_utils.get_stats_sig(posthoc_p[i]),
                    ha="center",
                )
                ax.axhline(
                    plot_ylim * (0.9 - i * 0.1),
                    xmin=sig_line_xs[i][0],
                    xmax=sig_line_xs[i][1],
                    color="black",
                    linestyle="-",
                    zorder=1,
                    lw=1,
                )
        ax.set_facecolor("whitesmoke")

    # Remove blanks
    if num_subplots < len(axes):
        for k in range(num_subplots, len(axes)):
            fig.delaxes(axes[k])
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
    plt.tight_layout()
    plt.show()
    plt.savefig(f"{out_dir}/{sfc_method}_{lm_lab}.svg", bbox_inches="tight")

# ---------------
# Plot difference
print(f"[Info] Plot INT-WMH", flush=True)
plot_ylim = 0.20
lm_lab = "INT-WMH"
fig, axes = plt.subplots(nrows=1, ncols=8, figsize=(8, 2.5))
axes = axes.flatten()
plt.subplots_adjust(wspace=0.2)
for i in range(num_subplots):
    cm_lab, cm_idx, _ = cm_arr[i]

    int_tbl = pd.read_csv(f"{out_dir}/{cm_lab}-cfc_{lm_labs[1]}.csv").assign(
        lm=lm_labs[1]
    )
    wmh_tbl = pd.read_csv(f"{out_dir}/{cm_lab}-cfc_{lm_labs[2]}.csv").assign(
        lm=lm_labs[2]
    )
    int_wmh_tbl = pd.concat([int_tbl, wmh_tbl])
    int_wmh_tbl = int_wmh_tbl[int_wmh_tbl["idx"] == cm_idx]

    glob_cmfc_tbl = int_wmh_tbl.pivot_table(
        values="rsq", index=["pid", "Group", "cm", "idx"], columns="lm"
    ).reset_index()
    glob_cmfc_tbl["rsq"] = glob_cmfc_tbl[lm_labs[1]] - glob_cmfc_tbl[lm_labs[2]]

    # ANOVA
    anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_ancova(
        glob_cmfc_tbl
    )

    comp_res.append(
        {
            "cm": cm_lab,
            "lm": "int-wmh",
            "anova_f": anova_f,
            "anova_p": anova_p,
            "p_group": p_group,
            "p_age": p_age,
            "posthoc_p_0": posthoc_p[0],
            "posthoc_p_1": posthoc_p[1],
            "posthoc_p_2": posthoc_p[2],
            mean_diff_lab[0]: meandiffs[0],
            mean_diff_lab[1]: meandiffs[1],
            mean_diff_lab[2]: meandiffs[2],
        }
    )
    ax = axes[i]
    ax = ana_utils.plot_ax_swarmplot(
        glob_cmfc_tbl, ax, i=0, x="Group", y="rsq", ylim=plot_ylim
    )

    ax.axhline(
        0.0,
        lw=1,
        color="black",
        linestyle="--",
        zorder=10,
        linewidth=1,
    )

    ax.set_ylabel(lm_lab)
    ax.set_ylim(-plot_ylim / 2, plot_ylim)
    ax.set_yticks(np.round(np.linspace(-plot_ylim / 2, plot_ylim, 4), 3))
    sns.despine(trim=True, ax=ax, offset=1)
    if i > 0:
        ax.spines["left"].set_visible(False)
        ax.set(yticks=[], ylabel="")

    ax.set_title(cm_lab[:2].upper())
    ax.set_xlabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if lm_idx == 2:
        ax.set_xticklabels(ana_utils.DISGROUPS_SHORT, rotation=90)
    else:
        ax.set_xticklabels([])
    sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
    for i in range(len(posthoc_p)):

        if p_group < 0.05 and posthoc_p[i] < 0.05:
            ax.text(
                0.5 + 0.5 * i,
                plot_ylim * (0.9 - i * 0.1),
                ana_utils.get_stats_sig(posthoc_p[i]),
                ha="center",
            )
            ax.axhline(
                plot_ylim * (0.9 - i * 0.1),
                xmin=sig_line_xs[i][0],
                xmax=sig_line_xs[i][1],
                color="black",
                linestyle="-",
                zorder=1,
                lw=1,
            )
    ax.set_facecolor("whitesmoke")
# Remove blanks
if num_subplots < len(axes):
    for k in range(num_subplots, len(axes)):
        fig.delaxes(axes[k])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
plt.tight_layout()
plt.show()
plt.savefig(f"{out_dir}/{sfc_method}_diff.svg", bbox_inches="tight")


comp_res_tbl = pd.DataFrame(comp_res)
comp_res_tbl.to_csv(out_dir + "/_comp_results.csv")

# -------------
# Draw heat map
# print("[Info] Draw heatmap for the meandiff")
# comp_res_tbl = pd.read_csv(out_dir + "/_comp_results.csv")
# for i in range(3):

#     comp_res_tbl[mean_diff_lab[i]] = np.where(
#         (comp_res_tbl["posthoc_p_" + str(i)] < 0.05) & (comp_res_tbl["p_group"] < 0.05),
#         comp_res_tbl[mean_diff_lab[i]],
#         0,
#     )

# for lm_idx, lm_lab in enumerate(lm_labs):
#     fig, axes = plt.subplots(figsize=(8, 1))
#     hm_dt = comp_res_tbl[comp_res_tbl.lm == lm_lab]
#     hm_dt = hm_dt[["cm", *mean_diff_lab]]
#     hm_dt = hm_dt.set_index("cm").transpose()
#     sns.heatmap(
#         cmap="RdBu_r",
#         vmax=0.05,
#         vmin=-0.05,
#         data=hm_dt,
#         linewidth=0.5,
#         cbar=False,
#     )
#     plt.xlabel("")

#     plt.show()
#     plt.savefig(f"{out_dir}/cms-sfc_anova_{lm_lab}.svg", bbox_inches="tight")
