import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import glob
import os

plt.style.use(['science', 'ieee'])

# データファイルのパス（ファイル名の並びに注意）
data_folder = '../../../../../../Results/4D-Var/Inverse/output/Ubend_X/average/velocity_profile_xy(z=bottom)/'
data_files = sorted(glob.glob(data_folder + 'velocity_phase*.csv'), key=lambda x: int(x.split('phase')[-1].split('.csv')[0]))
data_ref_files = sorted(glob.glob(data_folder + 'data_phase*.csv'), key=lambda x: int(x.split('phase')[-1].split('.csv')[0]))


# 保存先フォルダ
save_folder = r'./Ubend_inlet/'
os.makedirs(save_folder, exist_ok=True)

for i, (vel_file, data_file) in enumerate(zip(data_files, data_ref_files)):
    vel_data = pd.read_csv(vel_file)
    ref_data = pd.read_csv(data_file)

    # プロット設定
    plt.figure(figsize=(10, 6))
    
    # 各フェーズのプロット
    plt.plot(vel_data.iloc[:, 0], vel_data.iloc[:, 3], label=f'Assimilated', linewidth=5, linestyle='-', color='black')
    plt.plot(ref_data.iloc[:, 0], ref_data.iloc[:, 9], label=f'MRI', linewidth=5, linestyle='--', color='gray')

    # 軸ラベルと設定
    plt.xlabel('$x$ [m]', fontsize=26)
    plt.ylabel('Velocity ($x$) [m/s]', fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.ylim(-0.1,0.1)
    plt.legend(fontsize=28, loc='upper right')

    # ax = plt.gca()  # 現在の軸を取得
    # ax.spines['top'].set_linewidth(2)
    # ax.spines['right'].set_linewidth(2)
    # ax.spines['bottom'].set_linewidth(2)
    # ax.spines['left'].set_linewidth(2)

    # 図を保存
    save_path = f'{save_folder}velocity_x_inlet_phase{i}.png'
    plt.savefig(save_path, dpi=300)
    print(f'Saved: {save_path}')

    # プロットをクリア
    plt.close()


for i, (vel_file, data_file) in enumerate(zip(data_files, data_ref_files)):
    vel_data = pd.read_csv(vel_file)
    ref_data = pd.read_csv(data_file)

    # プロット設定
    plt.figure(figsize=(10, 6))
    
    # 各フェーズのプロット
    plt.plot(vel_data.iloc[:, 0], vel_data.iloc[:, 4], label=f'Assimilated', linewidth=4, linestyle='-', color='black')
    plt.plot(ref_data.iloc[:, 0], ref_data.iloc[:, 10], label=f'MRI', linewidth=4, linestyle='--', color='gray')

    # 軸ラベルと設定
    plt.xlabel('$x$ [m]', fontsize=26)
    plt.ylabel('Velocity ($y$) [m/s]', fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.ylim(0.0,0.7)
    plt.legend(fontsize=28, loc='upper right')

    # ax = plt.gca()  # 現在の軸を取得
    # ax.spines['top'].set_linewidth(2)
    # ax.spines['right'].set_linewidth(2)
    # ax.spines['bottom'].set_linewidth(2)
    # ax.spines['left'].set_linewidth(2)

    # 図を保存
    save_path = f'{save_folder}velocity_y_inlet_phase{i}.png'
    plt.savefig(save_path, dpi=300)
    print(f'Saved: {save_path}')

    # プロットをクリア
    plt.close()

for i, (vel_file, data_file) in enumerate(zip(data_files, data_ref_files)):
    vel_data = pd.read_csv(vel_file)
    ref_data = pd.read_csv(data_file)

    # プロット設定
    plt.figure(figsize=(10, 6))
    
    # 各フェーズのプロット
    plt.plot(vel_data.iloc[:, 0], vel_data.iloc[:, 5], label=f'Assimilated', linewidth=5, linestyle='-', color='black')
    plt.plot(ref_data.iloc[:, 0], ref_data.iloc[:, 11], label=f'MRI', linewidth=5, linestyle='--', color='gray')

    # 軸ラベルと設定
    plt.xlabel('$x$ [m]', fontsize=26)
    plt.ylabel('Velocity ($z$) [m/s]', fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.ylim(-0.1,0.1)
    plt.legend(fontsize=28, loc='upper right')

    # ax = plt.gca()  # 現在の軸を取得
    # ax.spines['top'].set_linewidth(2)
    # ax.spines['right'].set_linewidth(2)
    # ax.spines['bottom'].set_linewidth(2)
    # ax.spines['left'].set_linewidth(2)

    # 図を保存
    save_path = f'{save_folder}velocity_z_inlet_phase{i}.png'
    plt.savefig(save_path, dpi=300)
    print(f'Saved: {save_path}')

    # プロットをクリア
    plt.close()

