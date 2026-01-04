import ana_utils
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
import nilearn as nil
from nilearn import plotting
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import matplotlib.pyplot as plt

fig_out_dir = f"{ana_utils.PROJ_HOME}/data/lms_fig"
tbl_out_dir = f"{ana_utils.PROJ_HOME}/results/lms_stat"


def fdr_correction(pvals, alpha=0.05):
    """
    Perform Benjamini-Hochberg FDR correction on a list of p-values.

    Parameters:
    pvals (list or np.array): List or array of p-values to correct.
    alpha (float): Significance level for FDR correction.

    Returns:
    np.array: Array of corrected p-values.
    """
    _, corrected_pvals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    return corrected_pvals

def main():
    raw_res = pd.read_csv(f"{ana_utils.PROJ_HOME}/results/lr_net_lm_cfc/_comp_results.csv")
    new_res = pd.DataFrame()
    for lm in ana_utils.LM_LAB:
        for cm in ana_utils.CM_LAB_ARR:
            sub_res = raw_res[(raw_res["lm"] == lm) & (raw_res["cm"] == cm)]
            anova_pvals = sub_res["p_group"].values
            print(f"Performing FDR correction for lm: {lm}, cm: {cm}, {anova_pvals}")

            corrected_significant = fdr_correction(anova_pvals, alpha=0.05)
            sub_res["fdr_p"] = corrected_significant
            new_res = pd.concat([new_res, sub_res], ignore_index=True)
    new_res.to_csv(f"{ana_utils.PROJ_HOME}/results/lr_net_lm_cfc/_comp_results_with_fdr.csv", index=False)
    
    # Lr model only
    raw_res = pd.read_csv(f"{ana_utils.PROJ_HOME}/results/lr_net_lm_cfc/_comp_results_lr.csv")
    new_res = pd.DataFrame()
    for lm in ana_utils.LM_LAB:
        sub_res = raw_res[(raw_res["lm"] == lm)]
        anova_pvals = sub_res["p_group"].values
        print(f"Performing FDR correction for lm: {lm}, cm: {cm}, {anova_pvals}")

        corrected_significant = fdr_correction(anova_pvals, alpha=0.05)
        sub_res["fdr_p"] = corrected_significant
        new_res = pd.concat([new_res, sub_res], ignore_index=True)
    new_res.to_csv(f"{ana_utils.PROJ_HOME}/results/lr_net_lm_cfc/_comp_results_lr_with_fdr.csv", index=False)
    
if __name__ == "__main__":
    main()     