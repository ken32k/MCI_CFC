import ana_utils
import numpy as np
import pandas as pd
import os, sys
from joblib import Parallel, delayed

import matplotlib.pyplot as plt

print(f"Run {__file__}", flush=True)
print("-" * 20, flush=True)

# Create the empty output dir
output_dir = f"{ana_utils.PROJ_HOME}/results/mat_stat"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.system(f"find {output_dir}/* -type f -delete")

mat_sum_data = []
lm_labs = ana_utils.LM_LAB
for sub_id in ana_utils.SUB_ID_LIST:
    fc_mat_path = f"{ana_utils.PROJ_HOME}/data/fcs/{sub_id}_fc-wb_{ana_utils.ATLAS}.npy"
    sum_mat = np.load(fc_mat_path)[0].sum()
    mat_sum_data.append(
        {
            "pid": sub_id,
            "Group": ana_utils.SUB_ID_DIAG_DICT[sub_id],
            "cm": "fc",
            "lm": "wb",
            "value": sum_mat,
        }
    )

    for lm_lab in lm_labs:
        cm_mat_path = f"{ana_utils.PROJ_HOME}/data/cm/sc-wei/{sub_id}_{lm_lab}_sc-wei_{ana_utils.ATLAS}.npy"
        sum_mat = np.load(cm_mat_path)[0].sum()
        mat_sum_data.append(
            {
                "pid": sub_id,
                "Group": ana_utils.SUB_ID_DIAG_DICT[sub_id],
                "cm": "sc-wei",
                "lm": lm_lab,
                "value": sum_mat,
            }
        )


mat_tbl = pd.DataFrame(mat_sum_data)
mat_tbl.to_csv(output_dir + "/sub_mat_sum.csv")


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
mat_plot_tbl = mat_tbl[mat_tbl.cm == "fc"]
ax = ana_utils.plot_ax_swarmplot(mat_plot_tbl, ax=axes, i= 0, y="value",
                                  ylim=mat_plot_tbl["value"].max()*1.5)
plt.tight_layout()
plt.savefig(
    f"{output_dir}/mat-fc.png",
    bbox_inches="tight",
)

# plot sc
for lm_lab in lm_labs:
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))
    mat_plot_tbl = mat_tbl[(mat_tbl.cm == "sc-wei")& (mat_tbl.lm == lm_lab)]
    ax = ana_utils.plot_ax_swarmplot(mat_plot_tbl, ax=axes, i= 0, y="value",
                                    ylim=mat_plot_tbl["value"].max()*1.5)
    plt.tight_layout()
    plt.savefig(
        f"{output_dir}/mat-sc-{lm_lab}.png",
        bbox_inches="tight",
    )
