#!/bin/bash

# 要处理的文件夹
input_dir="/public_bme/data/v-baishw/CSVD/data/cm"
# 列出所有子文件夹
sub_dirs=("co-wei" "fg-wei" "mfpt-wei" "nav-pl" "pl-wei" "pt-wei" "sc-wei" "si-wei")
# 要搜索的文件模式（示例中包括了所有类型的文件）
file_patterns=("wb" "wmh" "int")

# 创建 CSV 文件并添加表头
echo -n "sub-" > result.csv
for dir in "${sub_dirs[@]}"
do
    echo -n ",$dir" >> result.csv
done
echo "" >> result.csv

# 遍历所有子文件夹和文件模式
for dir in "${sub_dirs[@]}"
do
    echo -n "$dir" >> result.csv
    for pattern in "${file_patterns[@]}"
    do
        count=$(find "$input_dir/$dir" -name "*_${pattern}_${dir}_schaefer-200.npy" | wc -l)
        echo -n ",$count" >> result.csv
    done
    echo "" >> result.csv
done