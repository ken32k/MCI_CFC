import ana_utils
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
import nilearn as nil
from nilearn import plotting
from scipy.stats import f_oneway
import statsmodels.api as sm
import matplotlib.pyplot as plt

fig_out_dir = f"{ana_utils.PROJ_HOME}/data/lms_fig"
tbl_out_dir = f"{ana_utils.PROJ_HOME}/results/lms_stat"

# main
for subject in ana_utils.SUB_ID_LIST:
    lm = np.load(f"{ana_utils.PROJ_HOME}/data/lms/{subject}_lm_schaefer-400.npy")
    fig, axes = plt.subplots(figsize=(8, 8))

    plotting.plot_matrix(mat=lm, 
                         colorbar=False, 
                         cmap="RdBu_r", 
                         vmin=-2, vmax=2)
    axes.axis("off")
    plt.show()
    plt.savefig(f"{fig_out_dir}/{subject}_lm_schaefer-400.png")

for lab in ana_utils.WMH_LAB_IDX[:-1]:
    node_stat = pd.DataFrame([])
    for subject in ana_utils.SUB_ID_LIST:
        # 读取患者的400的矩阵
        lm = np.load(f"{ana_utils.PROJ_HOME}/data/lms/{subject}_lm_schaefer-400.npy")

        # 按列累加，变成一个1*400的矩阵
        lm_count =  (lm == lab).sum(axis=0)

        # 创建临时 DataFrame 存储数据
        temp_df = pd.DataFrame(lm_count.reshape(1, -1))

        # 增加患者编号和分组信息
        temp_df['subject'] = subject
        temp_df['Group'] = ana_utils.SUB_ID_DIAG_DICT[subject]

        # 将临时 DataFrame 合并到 node_stat
        node_stat = pd.concat([node_stat, temp_df], ignore_index=True)
        
    node_stat.columns = [*['col' + str(i) for i in range(1, 401)], 'subject', 'Group']

    node_stat.to_csv(f"{tbl_out_dir}/subject_{ana_utils.WMH_LAB_VAL[lab]}_raw.csv")

    # node_stat = node_stat[~node_stat.subject]
    node_stat = node_stat.set_index(["Group","subject"])
    grouped_mean = node_stat.groupby('Group').mean()
    grouped_mean.to_csv(f"{tbl_out_dir}/subject_{ana_utils.WMH_LAB_VAL[lab]}_stat.csv")
