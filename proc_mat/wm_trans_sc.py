import os
import numpy as np

# 指定文件夹路径
folder_path = '/public_bme/data/v-baishw/CSVD/data/cm/sc-wei/'

# 获取文件夹下所有.npy文件
file_list = [f for f in os.listdir(folder_path) if f.endswith('.npy')]

# 遍历每个.npy文件
for file_name in file_list:
    try:
        file_path = os.path.join(folder_path, file_name)
        
        # 读取.npy文件
        data = np.load(file_path)
        
        # 修改shape为(1, 200, 200)
        data = data.reshape(1, 200, 200)
        
        # 保存修改后的数据，覆盖原文件
        np.save(file_path, data)
        
        print(f'{file_name} 处理完成',flush=True)
    except:
        print(file_name,flush=True)


print('所有文件处理完成')
