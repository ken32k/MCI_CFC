'''
Author: Ken32g ken32k@163.com
Date: 2024-11-26 22:30:54
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2024-11-26 22:30:55
FilePath: /csvd-sfc/analysis/plot_net_cfc.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import ana_utils
import numpy as np
import os
from nilearn import plotting
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

print("Plot regional cm function coupling.", flush=True)

for ana_method in ["corr","lr"]:
    csv_dir = f"{ana_utils.PROJ_HOME}/results/{ana_method}_net_lm_cfc"

    comp_res, mean_rsqs = [], []
    mean_diff_lab = ana_utils.mean_diff_lab
    disease_groups = ana_utils.DISGROUPS

    # Main loop
    for lm_lab in ana_utils.LM_LAB:
        for cm in ana_utils.CM_ARR[1:]:
            cm_lab, cm_idx, _ = cm
            reg_cfc_tbl = pd.read_csv(f"{csv_dir}/{cm_lab}-cfc_{lm_lab}.csv")
            uncorr_pval, sub_comp_res = [],  []
            for yeo_net_name in ana_utils.YEO7_DICT.keys():
                reg_cfc_net = reg_cfc_tbl[reg_cfc_tbl.net == yeo_net_name]
                mean_rsq = reg_cfc_net.groupby("Group")["rsq"].mean()
                mean_rsq_dict = {"cm": cm_lab, "lm": lm_lab, "net": yeo_net_name}
                for group in disease_groups:
                    mean_rsq_dict[group] = mean_rsq[group]
                mean_rsqs.append(mean_rsq_dict)
                anova_f, anova_p, p_group, p_age, posthoc_p, meandiffs = (
                    ana_utils.pg_anova(reg_cfc_net)
                )
                uncorr_pval.append(p_group)
                sub_comp_res.append(
                    {
                        "cm": cm_lab,
                        "lm": lm_lab,
                        "net": yeo_net_name,
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
            _, fdr_pval = fdrcorrection(uncorr_pval)
            for idx, pval in enumerate(fdr_pval):
                sub_comp_res[idx]["p_corr"] = fdr_pval[idx]
                comp_res.append(sub_comp_res[idx])


    comp_res_tbl = pd.DataFrame(comp_res)
    comp_res_tbl.to_csv(csv_dir + "/_comp_results.csv")
    meanrsq_tbl = pd.DataFrame(mean_rsqs)
    meanrsq_tbl.to_csv(csv_dir + "/_mean_rsq_results.csv")
