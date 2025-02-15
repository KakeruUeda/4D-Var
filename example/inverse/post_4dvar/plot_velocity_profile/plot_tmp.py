import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import os

plt.style.use(['science', 'ieee'])

# データファイルのパス
data_folder = '../../../../../../Results/4D-Var/Inverse/output/Ubend_X/bend_section1/'

da_file = 'tmp_DA.csv'
mri_file = 'tmp_MRI.csv'

# 保存先フォルダ
save_folder = './tmp/'
os.makedirs(save_folder, exist_ok=True)

# CSVデータの読み込み
da_data = pd.read_csv(da_file)
mri_data = pd.read_csv(mri_file)

# プロット設定
plt.figure(figsize=(10, 6))

# データプロット
plt.plot(da_data.iloc[:, 0], da_data.iloc[:, 4], label='Direct ($\mathbf{V}^{\mathrm{p}}$)', linewidth=5, linestyle='-', color='black')
plt.plot(mri_data.iloc[:, 3], mri_data.iloc[:, 1], label='MRI', linewidth=5, linestyle='--', color='gray')

# 軸ラベルと設定
plt.xlabel('$x$ [m]', fontsize=26)
plt.ylabel('Velocity (y) [m/s]', fontsize=26)
plt.xticks(fontsize=26)
plt.yticks(fontsize=26)
plt.legend(fontsize=28, loc='upper right')

# 図を保存
save_path = os.path.join(save_folder, 'velocity_comparison.png')
plt.savefig(save_path, dpi=300)
print(f'Saved: {save_path}')

# プロットをクリア
plt.close()