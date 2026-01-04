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

sfc_method = "corr_net"
output_dir = f"{ana_utils.PROJ_HOME}/results/{sfc_method}_lm_cfc"


mean_diff_lab = ana_utils.mean_diff_lab
cm_arr = ana_utils.CM_ARR[1:]
lm_labs = ana_utils.LM_LAB
net_lab_abbrs = ana_utils.net_abbr
net_labs = list(net_lab_abbrs.keys())
net_labs = [net_labs[0], net_labs[3], net_labs[6]]
comp_res = []
num_columns = len(cm_arr)
num_rows = 1

for net_idx, net_lab in enumerate(net_labs):
    # Plot all cfcs
    plot_ylim = 0.3

    for lm_idx, lm_lab in enumerate(lm_labs):
        print(f"[Info] Plot {lm_lab}", flush=True)
        fig, axes = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(3, 2.5))
        axes = axes.flatten()
        plt.subplots_adjust(left=0, wspace=0.1)
        for subplot_idx in range(num_columns):
            cm_lab, cm_idx, _ = cm_arr[subplot_idx]

            glob_cmfc_tbl = pd.read_csv(f"{output_dir}/{cm_lab}-cfc_{lm_lab}.csv")
            glob_cmfc_tbl = glob_cmfc_tbl[glob_cmfc_tbl.net == net_lab]
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
                    "net": net_lab,
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
            ax = axes[subplot_idx]
            ax = ana_utils.plot_ax_swarmplot(
                glob_cmfc_tbl,
                ax,
                i=0,
                x="Group",
                y="rsq",
                ylim=(0, plot_ylim),
                boxwidth=0.4,
                swarmsize=2,
            )
            # Set title
            if lm_idx == 0:
                ax.set_title(cm_lab[:2].upper())
            else:
                ax.set_title("")
            # spine
            ax.spines["top"].set_visible(False)
            sns.despine(trim=False, ax=ax, offset=1)

            # y axis
            ax.set_ylabel(["All", "INT", "WMH"][lm_idx])
            ax.set_ylim(0, plot_ylim)
            ax.set_yticks(np.round(np.linspace(0, plot_ylim, 4), 2))
            ax.spines["left"].set_visible(True)
            # plt.gca().yaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
            if net_idx >= 10 or subplot_idx > 0:
                ax.spines["left"].set_visible(False)
                ax.set(yticks=[], ylabel="")

            # x axis
            # ax.xaxis.set_visible(False)
            # ax.set_xlabel("")
            if lm_idx==2:
                ax.set_xticklabels(ana_utils.DISGROUPS_SHORT, rotation=90)
            else:
                ax.set_xticklabels([])

            # draw significance
            sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
            for p_idx in range(len(posthoc_p)):
                if p_group < 0.05 and posthoc_p[p_idx] < 0.05:
                    ax.text(
                        0.5 + 0.5 * p_idx,
                        plot_ylim * (0.9 - p_idx * 0.1),
                        ana_utils.get_stats_sig(posthoc_p[p_idx]),
                        ha="center",
                    )
                    ax.axhline(
                        plot_ylim * (0.9 - p_idx * 0.1),
                        xmin=sig_line_xs[p_idx][0],
                        xmax=sig_line_xs[p_idx][1],
                        color="black",
                        linestyle="-",
                        zorder=1,
                        lw=1,
                    )

        # plt.tight_layout()
        plt.savefig(
            f"{output_dir}/{sfc_method}_{net_lab}_{lm_lab}.png",
            dpi=300,
            bbox_inches="tight",
        )

    # ---------------
    # Plot difference
    print(f"[Info] Plot INT-WMH", flush=True)
    plot_ylim = 0.45
    lm_lab = "INT-WMH"
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(3, 2.5))
    axes = axes.flatten()
    plt.subplots_adjust(left=0, wspace=0.1)
    for subplot_idx in range(num_columns):
        cm_lab, cm_idx, _ = cm_arr[subplot_idx]

        int_tbl = pd.read_csv(f"{output_dir}/{cm_lab}-cfc_{lm_labs[1]}.csv").assign(
            lm=lm_labs[1]
        )
        wmh_tbl = pd.read_csv(f"{output_dir}/{cm_lab}-cfc_{lm_labs[2]}.csv").assign(
            lm=lm_labs[2]
        )
        int_wmh_tbl = pd.concat([int_tbl, wmh_tbl])
        int_wmh_tbl = int_wmh_tbl[int_wmh_tbl["idx"] == cm_idx]
        int_wmh_tbl = int_wmh_tbl[int_wmh_tbl.net == net_lab]

        glob_cmfc_tbl = int_wmh_tbl.pivot_table(
            values="rsq", index=["pid", "Group", "cm", "idx"], columns="lm"
        ).reset_index()
        glob_cmfc_tbl["rsq"] = glob_cmfc_tbl[lm_labs[1]] - glob_cmfc_tbl[lm_labs[2]]

        # ANOVA
        anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_anova(
            glob_cmfc_tbl
        )

        comp_res.append(
            {
                "cm": cm_lab,
                "lm": "int-wmh",
                "net": net_lab,
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
        ax = axes[subplot_idx]
        ax = ana_utils.plot_ax_swarmplot(
            glob_cmfc_tbl,
            ax,
            i=0,
            x="Group",
            y="rsq",
            ylim=plot_ylim,
            boxwidth=0.4,
            swarmsize=2,
        )
        # draw zero lines
        ax.axhline(
            0.0,
            lw=1,
            color="black",
            linestyle="--",
            zorder=10,
            linewidth=1,
        )

        ax.set_title("")

        # spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # y axis
        ax.set_ylabel(lm_lab)
        ax.set_ylim(-plot_ylim / 2, plot_ylim)

        ax.set_yticks(np.round(np.linspace(-plot_ylim / 2, plot_ylim, 4), 3))
        sns.despine(trim=True, ax=ax, offset=1)
        if net_idx >= 10 or subplot_idx > 0:
            ax.spines["left"].set_visible(False)
            ax.set(yticks=[], ylabel="")

        ax.set_xlabel("")
        ax.set_xticklabels(ana_utils.DISGROUPS_SHORT, rotation=90)

        # draw significance lines
        sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
        for p_idx in range(len(posthoc_p)):

            if p_group < 0.05 and posthoc_p[p_idx] < 0.05:
                ax.text(
                    0.5 + 0.5 * p_idx,
                    plot_ylim * 1.5 * (0.9 - p_idx * 0.1) - 0.5 * plot_ylim,
                    ana_utils.get_stats_sig(posthoc_p[p_idx]),
                    ha="center",
                )
                ax.axhline(
                    plot_ylim * 1.5 * (0.9 - p_idx * 0.1) - 0.5 * plot_ylim,
                    xmin=sig_line_xs[p_idx][0],
                    xmax=sig_line_xs[p_idx][1],
                    color="black",
                    linestyle="-",
                    zorder=1,
                    lw=1,
                )

    plt.gca().yaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
    # plt.tight_layout()
    plt.savefig(
        f"{output_dir}/{sfc_method}_{net_lab}_diff.png", dpi=300, bbox_inches="tight"
    )


comp_res_tbl = pd.DataFrame(comp_res)
comp_res_tbl.to_csv(output_dir + "/_comp_results.csv")

# lm_labs = [*lm_labs, "diff"]
# Combine subplots
grid_figures = []
for lm_lab in lm_labs:
    for net_lab in net_labs:
        grid_figures.append(f"{output_dir}/{sfc_method}_{net_lab}_{lm_lab}.png")

ana_utils.plot_create_image_grid(
    grid_figures,
    len(lm_labs),
    len(net_labs),
    30,
    f"{output_dir}/grid_rada.png",
)
