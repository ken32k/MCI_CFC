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
from joblib import Parallel, delayed


def par_net_cfc_lrm(cm_list, sub_id, lm_lab):
    """Parallel calculate corr with global lrm"""
    sub_lrm_result = []
    for yeo_net_name in ana_utils.YEO7_DICT.keys():
        # Get communication model matrices
        cm_mats = []
        for cm in cm_list:
            cm_lab, cm_idx = cm[0], cm[1]
            cm_mat = np.load(
                f"{ana_utils.PROJ_HOME}/data/cm/{cm_lab}/{sub_id}_{lm_lab}_{cm_lab}_{ana_utils.ATLAS}.npy"
            )
            cm = cm_mat[cm_idx][
                np.ix_(ana_utils.YEO7_DICT[yeo_net_name],
                       ana_utils.YEO7_DICT["Default"])]
            cm_mats.append(cm)
        # Get fc
        fc_mat = np.load(
            f"{ana_utils.PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
        )[np.ix_(ana_utils.YEO7_DICT[yeo_net_name],
                ana_utils.YEO7_DICT["Default"])]

        # Get rsq, pval and dominance
        try:
            rsq_val, p_val, dominance, _ = ana_utils.sub_mat_linear_reg(
                np.power(fc_mat,2), np.array(cm_mats)
            )
        except Exception as e:
            rsq_val, p_val = np.nan, np.nan
            dominance = [np.nan for i in range(5)]

        sub_net_res = {
            "pid": sub_id,
            "Group": ana_utils.SUB_ID_DIAG_DICT[sub_id],
            "cm": "lrm",
            "idx": 0,
            "rsq": rsq_val,
            "pval": p_val,
            "net": yeo_net_name,
        }
        # Add dominance
        for i, cm in enumerate(cm_list):
            cm_lab, cm_idx = cm[0], cm[1]
            sub_net_res[cm_lab] = dominance[i]
        sub_lrm_result.append(sub_net_res)
    return sub_lrm_result


def get_cfc_lrm(cm_list, lm_lab):
    """Perform correlation and save csv"""
    # Get subject list
    sub_cm_list = pd.read_csv(f"{ana_utils.PROJ_HOME}/data/sub_list.csv")[
        "pid"
    ].to_list()

    # Main parallel loop
    cfc_results =np.array(
        Parallel(n_jobs=-1)(
            delayed(par_net_cfc_lrm)(cm_list, sub_id, lm_lab) for sub_id in sub_cm_list
        )
    )

    return cfc_results


print(f"Run {__file__}", flush=True)
print("-" * 20, flush=True)

# Create the output dir
output_dir = f"{ana_utils.PROJ_HOME}/results/lr_net_lm_cfc"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# # Delete existing files
os.system(f"find {output_dir}/* -type f -delete")

mean_diff_lab = ana_utils.mean_diff_lab
cm_list = ana_utils.CM_ARR[1:]

# Get RSQ with Linear regression model
for i, lm_lab in enumerate(ana_utils.LM_LAB):
    cfc_results = get_cfc_lrm(cm_list, lm_lab)
    
    flattened_list = [item for sublist in cfc_results for item in sublist]
    cfc_results_tbl = pd.json_normalize(flattened_list)
    cfc_results_tbl.to_csv(
        f"{output_dir}/lrcms-cfc_{lm_lab}.csv",
        index=False,
    )

    for cm in cm_list:
        cm_lab, cm_idx = cm[0], cm[1]
        cfc_results_tbl.loc[:, "rsq"] = cfc_results_tbl[cm_lab]
        cm_dom_tbl = cfc_results_tbl[["pid", "Group", "cm", "idx", "rsq", "net"]]
        cm_dom_tbl.loc[:, "idx"] = cm_idx
        cm_dom_tbl.to_csv(
            f"{output_dir}/{cm_lab}-cfc_{lm_lab}.csv",
            index=False,
        )


for yeo_net_name in ana_utils.YEO7_DICT.keys():
    # Anova
    comp_res = []
    # Plot
    num_columns = len(ana_utils.LM_LAB)  # 分成的列数
    fig, axes = plt.subplots(nrows=1, ncols=num_columns, figsize=(8, 6))
    axes = axes.flatten()

    for i, lm_lab in enumerate(ana_utils.LM_LAB):
        # plot
        ax = axes[i]
        lm_res = pd.read_csv(f"{output_dir}/lrcms-cfc_{lm_lab}.csv")
        lm_res = lm_res[lm_res.net == yeo_net_name]
        anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = ana_utils.pg_anova(
            lm_res
        )
        comp_res.append(
            {
                "cm": 0,
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
        ax = ana_utils.plot_ax_swarmplot(lm_res, ax, i)
        ax.set_xlabel("")
        ax.set_ylabel(lm_lab)
        mean_rsq = lm_res.groupby("Group")["rsq"].mean()
        m = 0
        # Draw mean line
        for group, mean in mean_rsq.items():
            ax.axhline(
                mean,
                xmin=m / 3 + 1 / 12,
                xmax=m / 3 + 1 / 12 + 1 / 6,
                color="black",
                linestyle="-",
                zorder=10,
            )  
            m += 1
        wmh_ylim = lm_res.rsq.max() * 1.6
        ax.set_yticks(np.round(np.linspace(0, wmh_ylim, 5), 3))
        if p_group < 0.05:
            ax.text(1, wmh_ylim * 0.9, ana_utils.get_stats_sig(anova_p), ha="center")
        ax.set_yticks(np.round(np.linspace(0, wmh_ylim, 5), 2))
    plt.tight_layout()
    plt.show()
    plt.savefig(
        f"{output_dir}/lmr_wb-sfc_{lm_lab}-{yeo_net_name}.png",
        bbox_inches="tight",
    )
import lr_net_lm_plot
os.system(f"mv {output_dir} {output_dir}-Default")