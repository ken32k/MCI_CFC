#!/usr/bin/env python3
import ana_utils
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp

import pingouin as pg
from matplotlib import colors as mcolors

print("Run ", __file__, flush=True)

output_dir = f"{ana_utils.PROJ_HOME}/results/lr_glob_lm_cfc"
print("Output dir:", output_dir, flush=True)

colors_hex = [mcolors.to_hex(color) for color in ana_utils.PAL4]
cm_arr_hex = [mcolors.to_hex(color) for color in ana_utils.CM_ARR_COLOR]
mean_diff_lab = ana_utils.mean_diff_lab
disease_groups = ana_utils.DISGROUPS
lm_labs = ana_utils.LM_LAB
CM_LAB_ARR = ana_utils.CM_LAB_ARR

# Merge cfc data for five LMs
merged_data = pd.concat(
    [
        pd.read_csv(f"{output_dir}/lrcms-cfc_{lm_lab}.csv").assign(lm=lm_lab)
        for lm_lab in lm_labs
    ],
    ignore_index=True,
)

# Calculate merged_diff (INT - WMH for each pid and Group)
int_subtbl = merged_data[merged_data.lm == lm_labs[1]].set_index(["pid", "Group"])
wmh_subtbl = merged_data[merged_data.lm == lm_labs[2]].set_index(["pid", "Group"])
merged_diff = (int_subtbl[CM_LAB_ARR] - wmh_subtbl[CM_LAB_ARR]).reset_index()

# Swarmplot for global R2 of three LMs
ylim = 0.28
fig, axes = plt.subplots(nrows=1, ncols=len(lm_labs), figsize=(2.5, 5))
axes = axes.flatten()
comp_res = []

for i, lm_lab in enumerate(lm_labs):
    ax = axes[i]
    lr_res_tbl = merged_data[merged_data.lm == lm_lab]
    anova_f, _, p_group, _, posthoc_p, meandiffs = ana_utils.pg_anova(lr_res_tbl)

    ana_utils.plot_ax_swarmplot(
        lr_res_tbl, ax, i, x="Group", y="rsq", ylim=0.24, swarmsize=3
    )
    ax.set_ylabel("")
    ax.set_yticks(np.round(np.linspace(0, ylim, 5), 2))
    ax.set_yticklabels(ax.get_yticks(), rotation=90)

    if i > 0:
        ax.spines["left"].set_color("None")
        ax.set(yticks=[], ylabel="")
    if ax.get_ylim()[0] > ax.get_ylim()[1]:
        ax.invert_yaxis()
    sns.despine(trim=True, ax=ax)
    sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
    for j in range(len(posthoc_p)):
        if p_group < 0.05 and posthoc_p[j] < 0.05:
            ax.text(
                0.5 + 0.5 * j,
                ylim * (0.9 - j * 0.1),
                ana_utils.get_stats_sig(posthoc_p[j]),
                ha="center",
            )
            ax.axhline(
                ylim * (0.9 - j * 0.1),
                xmin=sig_line_xs[j][0],
                xmax=sig_line_xs[j][1],
                color="black",
                linestyle="-",
                lw=1,
                zorder=1,
            )
    comp_res.append(
        {
            "cm": "lr_glo",
            "lm": lm_lab,
            "anova_f": anova_f,
            "anova_p": 0,
            "p_group": p_group,
            "p_age": 0,
            "posthoc_p_0": posthoc_p[0],
            "posthoc_p_1": posthoc_p[1],
            "posthoc_p_2": posthoc_p[2],
            mean_diff_lab[0]: meandiffs[0],
            mean_diff_lab[1]: meandiffs[1],
            mean_diff_lab[2]: meandiffs[2],
        }
    )
pd.DataFrame(comp_res).to_csv(f"{output_dir}/_comp_results_lr.csv")

plt.tight_layout()
plt.savefig(f"{output_dir}/lr_glo.svg", bbox_inches="tight", dpi=300)
# 


# ANOVA for each CM and LM
comp_res = []
for lm_lab in lm_labs:
    for cm_lab in CM_LAB_ARR:
        lr_glob_cmfc_tbl = merged_data[merged_data.lm == lm_lab]
        anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_anova(
            lr_glob_cmfc_tbl, cm_lab
        )
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
pd.DataFrame(comp_res).to_csv(f"{output_dir}/_comp_results.csv")

# Plot radar charts
def plot_radar(ax, values, angles, color, label, linestyle="solid"):
    ax.plot(angles, values, linewidth=1, linestyle=linestyle, color=color, label=label)
    ax.set_facecolor("white")

def get_angles(N):
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    return angles + angles[:1]

def get_rsq_mean(tbl, groupby="Group", index="pid"):
    return tbl.set_index(index).groupby(groupby).mean()

def plot_radar_for_tbl(tbl, title, ylim, colors, labels, angles, ax, linestyle="solid"):
    for idx, label in enumerate(labels):
        values = tbl.loc[label].values
        values = np.concatenate((values, [values[0]]))
        plot_radar(ax, values, angles, colors[idx], label, linestyle)
    ax.set_title(title)
    ax.set_yticks(np.round(np.linspace(0, ylim, 3), 3))
    plt.xticks(angles[:-1], ana_utils.CM_ARR_SHORT)
    locs, xticklabels = plt.xticks()
    [plt.setp(xticklabels[i], color=cm_arr_hex[i]) for i in range(len(xticklabels))]

# All fibers
ylim = 0.03
lm_lab = lm_labs[0]
lr_res_tbl = merged_data[merged_data.lm == lm_lab][["Group", "pid", *[item[0] for item in ana_utils.CM_ARR[1:]]]]
lr_rsq_mean = get_rsq_mean(lr_res_tbl)
categories = list(lr_rsq_mean.columns)
N = len(categories)
angles = get_angles(N)
fig, ax = plt.subplots(figsize=(3, 3), subplot_kw={"polar": True})
plot_radar_for_tbl(lr_rsq_mean, "All", ylim, colors_hex, disease_groups, angles, ax)
plt.savefig(f"{output_dir}/rada_lm-{lm_lab}.svg", bbox_inches="tight")

# INT & WMH
fig, ax = plt.subplots(figsize=(3, 3), subplot_kw={"polar": True})
for lm_lab, linestyle in zip(lm_labs[1:3], ["solid", ":"]):
    lr_res_tbl = merged_data[merged_data.lm == lm_lab][["Group", *[item[0] for item in ana_utils.CM_ARR[1:]]]]
    lr_rsq_mean = lr_res_tbl.groupby("Group").mean()
    plot_radar_for_tbl(lr_rsq_mean, "", ylim, colors_hex, disease_groups, get_angles(len(lr_rsq_mean.columns)), ax, linestyle)
ax.set_title("INT & WMH")
plt.savefig(f"{output_dir}/rada_lm-dec.svg", bbox_inches="tight")

# Difference
lr_res_tbl = merged_diff[["Group", "pid", *[item[0] for item in ana_utils.CM_ARR[1:]]]]
lr_rsq_mean = get_rsq_mean(lr_res_tbl)
categories = list(lr_rsq_mean.columns)
N = len(categories)
angles = get_angles(N)
fig, ax = plt.subplots(figsize=(3, 3), subplot_kw={"polar": True})
plot_radar_for_tbl(lr_rsq_mean, "INT-WMH", 0.02, colors_hex, disease_groups, angles, ax)
ax.set_yticks(np.round(np.linspace(-0.01, 0.02, 3), 2))
plt.savefig(f"{output_dir}/rada_lm-intwmh.svg", bbox_inches="tight")

# Pearson correlation and ANOVA
pearson_stat, comp_res = [], []
pearson_mat = np.zeros((len(CM_LAB_ARR), len(CM_LAB_ARR)))
for i, cm1 in enumerate(CM_LAB_ARR):
    for j, cm2 in enumerate(CM_LAB_ARR):
        for lm_lab in lm_labs:
            wb_tbl = merged_data[merged_data.lm == lm_lab]
            r, p = sp.stats.pearsonr(wb_tbl[cm1], wb_tbl[cm2])
            pearson_stat.append({"cat": lm_lab, "int": cm1, "wmh": cm2, "r": r, "p": p})

        int_wmh_tbl = merged_data[merged_data.lm.isin(lm_labs[1:])]
        res_pt = pg.pairwise_ttests(int_wmh_tbl, dv=cm1, between="Group", subject="pid", within="lm")
        print(res_pt)

        int_wmh_tbl = int_wmh_tbl.set_index("pid").pivot_table(index=["pid", "Group"], values=CM_LAB_ARR, columns="lm")
        fig, ax = plt.subplots(figsize=(3, 3))
        sns.regplot(x=int_wmh_tbl[cm1]["int"], y=int_wmh_tbl[cm2]["wmh"], ax=ax)
        r, p = sp.stats.pearsonr(int_wmh_tbl[cm1]["int"], int_wmh_tbl[cm2]["wmh"])
        pearson_stat.append({"cat": "int-wmh", "int": cm1, "wmh": cm2, "r": r, "p": p})
        if p < 0.05:
            pearson_mat[i][j] = r
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{cm1}-{cm2}_reg.svg", bbox_inches="tight")

        int_wmh_tbl["cm1-2"] = int_wmh_tbl[cm1]["int"] - int_wmh_tbl[cm2]["wmh"]
        fig, ax = plt.subplots(figsize=(2, 4))
        print(ana_utils.pg_anova(int_wmh_tbl.reset_index(), "cm1-2")[2])
        ana_utils.plot_ax_swarmplot(int_wmh_tbl, swarmsize=3, ax=ax, x="Group", y="cm1-2")
        anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_anova(
            int_wmh_tbl.reset_index(), "cm1-2"
        )
        comp_res.append(
            {
                "cm1": cm1,
                "cm2": cm2,
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
        ax.set_yticks(np.round(np.linspace(0, 2, 3), 3))
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{cm1}-{cm2}.svg", bbox_inches="tight")

# Draw heatmap
fig, ax = plt.subplots(figsize=(4, 4))
sns.heatmap(
    pearson_mat,
    cmap="RdBu_r",
    vmax=1,
    vmin=-1,
    square=True,
    xticklabels=[cm[:2].upper() for cm in CM_LAB_ARR],
    yticklabels=[cm[:2].upper() for cm in CM_LAB_ARR],
)
plt.tight_layout()
plt.savefig(f"{output_dir}/corr_heat.svg", bbox_inches="tight")

pd.DataFrame(pearson_stat).to_csv(f"{output_dir}/pearson_corr.csv")
pd.DataFrame(comp_res).to_csv(f"{output_dir}/int_wmh_comp.csv")
