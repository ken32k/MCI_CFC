import os
import numpy as np
import pandas as pd
import ana_utils
import matplotlib.pyplot as plt
import pingouin as pg
import matplotlib.colors as mcolors
import seaborn as sns
from PIL import Image, ImageOps

print("Run ", __file__, flush=True)

output_dir = f"{ana_utils.PROJ_HOME}/results/lr_net_lm_cfc"

dis_groups=ana_utils.DISGROUPS
colors_hex = [mcolors.to_hex(color) for color in ana_utils.PAL4]
cm_arr_hex = [mcolors.to_hex(color) for color in ana_utils.CM_ARR_COLOR]
cm_lab_labs = ana_utils.CM_LAB_ARR
mean_diff_lab = ana_utils.mean_diff_lab
net_lab_abbrs = ana_utils.net_abbr
net_labs = net_lab_abbrs.keys()
lm_labs = ana_utils.LM_LAB
data_list, file_list = [], []

comp_res = []

plt.rcParams.update(
    {"font.size": 15}
)  # Change the font size (e.g., 12) as per your requirement
# Merge all csv files
merged_data = pd.DataFrame()
for i, lm_lab in enumerate(lm_labs):
    data = pd.read_csv(f"{output_dir}/lrcms-cfc_{lm_lab}.csv")
    data = data.assign(lm=lm_lab)
    merged_data = pd.concat(
        [merged_data, data],
        ignore_index=True,
    )
    
merged_data
# Calculate merged_diff
merged_diff = pd.DataFrame()
for net_lab in net_labs:
    int_subtbl = merged_data[merged_data.lm == lm_labs[1]].set_index(
        ["pid", "net", "Group"]
    )
    wmh_subtbl = merged_data[merged_data.lm == lm_labs[2]].set_index(
        ["pid", "net", "Group"]
    )

    df_tbl = int_subtbl[cm_lab_labs] - wmh_subtbl[cm_lab_labs]
    merged_diff = pd.concat(
        [merged_diff, df_tbl.reset_index()],
        ignore_index=True,
    )

# Plot LR rsq
for lm_idx, lm_lab in enumerate(lm_labs):
    fig, axes = plt.subplots(nrows=1, ncols=len(net_labs), figsize=(8, 3))
    plt.subplots_adjust(wspace=0.02)
    axes = axes.flatten()
    for net_idx, net_lab in enumerate(net_labs):
        # plot
        ax = axes[net_idx]
        lr_res_tbl = merged_data[
            (merged_data.lm == lm_lab) & (merged_data.net == net_lab)
        ]
        anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_anova(
            lr_res_tbl
        )
        comp_res.append(
            {
                "cm": "wholelr",
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
      
        ylim = 0.5
        ax = ana_utils.plot_ax_swarmplot(lr_res_tbl, ax, net_idx, x="Group", y="rsq", swarmsize=3)
        ax.set_ylabel(["All", "INT", "WMH"][lm_idx])

        ax.set_title(net_lab_abbrs[net_lab])

        ax.set_yticks(np.round(np.linspace(0, ylim, 5), 2))
        if net_idx > 0:
            ax.spines["left"].set_color("None")
            ax.set(yticks=[], ylabel="")
        if ax.get_ylim()[0] > ax.get_ylim()[1]:
            ax.invert_yaxis()
        sns.despine(trim=True, ax=ax)

        sig_line_xs = [(1 / 6, 1 / 2), (1 / 6, 5 / 6), (1 / 2, 5 / 6)]
        for i in range(len(posthoc_p)):
            if p_group < 0.05 and posthoc_p[i] < 0.05:
                ax.text(
                    0.5 + 0.5 * i,
                    ylim * (0.9 - i * 0.1),
                    ana_utils.get_stats_sig(posthoc_p[i]),
                    ha="center",
                )
                ax.axhline(
                    ylim * (0.9 - i * 0.1),
                    xmin=sig_line_xs[i][0],
                    xmax=sig_line_xs[i][1],
                    color="black",
                    linestyle="-",
                    zorder=1,
                )
    plt.tight_layout()
    plt.savefig(f"{output_dir}/lr_net-{lm_lab}.svg", bbox_inches="tight")

comp_res_tbl = pd.DataFrame(comp_res)
comp_res_tbl.to_csv(output_dir + "/_comp_results_lr.csv")




# =========
# Plot rada
print("")
rada_fig_size = (3, 3)
ylim = 0.04
for net_idx, net_lab in enumerate(net_labs):

    # Whole brain
    lm_lab = lm_labs[0]
    lr_res_tbl = merged_data[(merged_data.lm == lm_lab) & (merged_data.net == net_lab)]
    lr_res_tbl = lr_res_tbl[["Group", "pid", *cm_lab_labs]]
    lr_rsq_mean = lr_res_tbl.set_index("pid").groupby("Group").mean()

    categories = list(lr_rsq_mean.columns)
    N = len(categories)
    fig, ax = plt.subplots(figsize=rada_fig_size, subplot_kw={"polar": True})
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    for idx, dis_group in enumerate(lr_rsq_mean.index):
        values = lr_rsq_mean.loc[dis_group].values
        values = np.concatenate((values, [values[0]]))

        ax.plot(
            angles,
            values,
            linewidth=1,
            linestyle="solid",
            color=colors_hex[dis_groups.index(dis_group)],
            label=dis_group,
        )
        ax.set_facecolor("white")
    pvals_list = []
    ax.set_yticks(np.round(np.linspace(0, ylim, 3), 3))

    plt.xticks(angles[:-1], ana_utils.CM_ARR_SHORT)
    locs, labels = plt.xticks()
    [plt.setp(labels[i], color=cm_arr_hex[i]) for i in range(len(labels))]
    plt.savefig(
        f"{output_dir}/rada_lm-{lm_lab}_net-{net_lab}.png",
        bbox_inches="tight",
        dpi=300,
    )

    # ------------------------
    # INT AND WMH
    # INT
    fig, ax = plt.subplots(figsize=rada_fig_size, subplot_kw={"polar": True})

    lm_lab = lm_labs[1]
    lr_res_tbl = merged_data[(merged_data.lm == lm_lab) & (merged_data.net == net_lab)]
    lr_res_tbl = lr_res_tbl[["Group", *cm_lab_labs]]
    lr_rsq_mean = lr_res_tbl.groupby("Group").mean()
    categories = list(lr_rsq_mean.columns)
    N = len(categories)

    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    for idx, dis_group in enumerate(lr_rsq_mean.index):
        values = lr_rsq_mean.loc[dis_group].values
        values = np.concatenate((values, [values[0]]))

        ax.plot(
            angles,
            values,
            linewidth=1,
            linestyle="solid",
            color=colors_hex[dis_groups.index(dis_group)],
            label=dis_group,
        )
    # --------------------------
    # WMH
    lm_lab = lm_labs[2]
    lr_res_tbl = merged_data[(merged_data.lm == lm_lab) & (merged_data.net == net_lab)]
    lr_res_tbl = lr_res_tbl[["Group", *cm_lab_labs]]
    lr_rsq_mean = lr_res_tbl.groupby("Group").mean()
    categories = list(lr_rsq_mean.columns)
    N = len(categories)

    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    for idx, dis_group in enumerate(lr_rsq_mean.index):
        values = lr_rsq_mean.loc[dis_group].values
        values = np.concatenate((values, [values[0]]))

        ax.plot(
            angles,
            values,
            linewidth=1,
            linestyle=":",
            color=colors_hex[dis_groups.index(dis_group)],
            label=dis_group,
        )
        ax.set_facecolor("white")
        ax.set_title("")

        ax.set_yticks(np.round(np.linspace(0, ylim, 3), 3))
    # ax.set_title("INT & WMH" if net_idx == 0 else "")
    plt.xticks(angles[:-1], ana_utils.CM_ARR_SHORT)
    locs, labels = plt.xticks()
    [plt.setp(labels[i], color=cm_arr_hex[i]) for i in range(len(labels))]
    plt.savefig(
        f"{output_dir}/rada_lm-dec_net-{net_lab}.png", bbox_inches="tight", dpi=300
    )
    # ---------------------
    # Difference

    lr_res_tbl = merged_diff[merged_diff.net == net_lab]
    lr_res_tbl = lr_res_tbl[["Group", "pid", *cm_lab_labs]]
    lr_rsq_mean = lr_res_tbl.set_index("pid").groupby("Group").mean()

    categories = list(lr_rsq_mean.columns)
    N = len(categories)
    fig, ax = plt.subplots(figsize=rada_fig_size, subplot_kw={"polar": True})
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    for idx, dis_group in enumerate(lr_rsq_mean.index):
        values = lr_rsq_mean.loc[dis_group].values
        values = np.concatenate((values, [values[0]]))

        ax.plot(
            angles,
            values,
            linewidth=1,
            linestyle="solid",
            color=colors_hex[dis_groups.index(dis_group)],
            label=dis_group,
        )
        ax.set_facecolor("white")
        ax.set_title("")
    # ax.set_title("INT-WMH" if net_idx == 0 else "")
    ax.set_yticks(np.round(np.linspace(-0.02, 0.02, 3), 2))
    plt.xticks(angles[:-1], ana_utils.CM_ARR_SHORT)
    locs, labels = plt.xticks()
    [plt.setp(labels[i], color=cm_arr_hex[i]) for i in range(len(labels))]
    plt.savefig(
        f"{output_dir}/rada_lm-intwmh_net-{net_lab}.png", bbox_inches="tight", dpi=300
    )


# Combine subplots
grid_figures = []

for lb in ["wb", "dec", "intwmh"]:
    for net in net_labs:
        grid_figures.append(f"{output_dir}/rada_lm-{lb}_net-{net}.png")
ana_utils.plot_create_image_grid(
    grid_figures,
    3,
    len(net_labs),
    20,
    f"{output_dir}/grid_rada.png",
    # [ana_utils.net_abbr[net] for net in net_labs],
    # [lm_labs]
)
