"""
Author: Ken32g ken32k@163.com
Date: 2024-11-27 12:58:24
LastEditors: Ken32g ken32k@163.com
LastEditTime: 2024-11-27 12:59:28
FilePath: /csvd-sfc/WMHseg/mci_wmh_summary_lst.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
"""

import pandas as pd
import numpy as np
import os, time

# Set the directory for the CSV files
HOME = os.getenv("CSVD_DIR")
bids_dir = f"{HOME}/csvd_mica_bids"
sub_list_beforpsm = f"{HOME}/data/sub_list_beforePSM.csv"
output_csv_path = f"{HOME}/data/sub_list_beforePSM_WMHstat.csv"

# Load the sub_list_beforePSM.csv into a DataFrame
sub_list_df = pd.read_csv(sub_list_beforpsm)

# Initialize an empty DataFrame to store WMH stats
wmh_stat = pd.DataFrame()

# Iterate over each subject directory
for subdir in os.listdir(bids_dir):
    sub_wmh_file = os.path.join(bids_dir, subdir, "wmh", "lesion_stats.csv")

    if os.path.exists(sub_wmh_file):
        # Load the lesion stats into a DataFrame
        sub_wmh_data = pd.read_csv(sub_wmh_file)

        # Add a new column for the subject ID
        sub_wmh_data["pid"] = subdir

        # Append the data to wmh_stat DataFrame
        wmh_stat = pd.concat([wmh_stat, sub_wmh_data], ignore_index=True)

# Merge the WMH stats with the subject list on "sub-id"
merged_data = pd.merge(sub_list_df, wmh_stat, on="pid", how="left")
merged_data = merged_data[merged_data["Lesion_Volume"] > 0]
merged_data["logwmh"] = np.log10(merged_data["Lesion_Volume"])
# Display the merged data
merged_data.to_csv(output_csv_path)
print(f"{output_csv_path}, created at {time.ctime(os.path.getctime(output_csv_path))}")
