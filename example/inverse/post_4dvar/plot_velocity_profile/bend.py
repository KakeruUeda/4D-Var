import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import glob
import os

plt.style.use(['science', 'ieee'])

# データファイルのパス（ファイル名の並びに注意）
data_folder = '../../../../../../Results/4D-Var/Inverse/output/Ubend_X/bend_section1/'
# ソート用の安全な関数を定義
def extract_phase_number(file_name):
    try:
        # phase番号を抽出
        return int(file_name.split('phase_')[-1].split('.csv')[0])
    except (IndexError, ValueError):
        print(f"Skipping invalid file: {file_name}")
        return float('inf')  # 無効なファイルは最後に配置

# velocity_y_inlet_phase_* ファイルの取得とソート
data_files = sorted(
    glob.glob(data_folder + 'velocity_y_inlet_phase_*.csv'),
    key=extract_phase_number
)

# data_y_inlet_phase_* ファイルの取得とソート
data_ref_files = sorted(
    glob.glob(data_folder + 'data_y_inlet_phase_*.csv'),
    key=extract_phase_number
)

direct_files = sorted(
    glob.glob(data_folder + 'direct_assimilated_velocity_y_inlet_phase_*.csv'),
    key=extract_phase_number
)

# 保存先フォルダ
save_folder = r'./Ubend_bend/'
os.makedirs(save_folder, exist_ok=True)

for i, (vel_file, data_file, direct_file) in enumerate(zip(data_files, data_ref_files, direct_files)):
    vel_data = pd.read_csv(vel_file)
    ref_data = pd.read_csv(data_file)
    direct_data = pd.read_csv(direct_file)

    # プロット設定
    plt.figure(figsize=(10, 6))
    
    # 各フェーズのプロ
    plt.plot(vel_data.iloc[:, 0], vel_data.iloc[:, 3], label=f'Assimilated', linewidth=5, linestyle='-', color='black')
    #plt.plot(vel_data.iloc[:, 0], direct_data.iloc[:, 3], label=f'Direct (Assimilated BC)', linewidth=5, linestyle='-', color='blue')
    plt.plot(ref_data.iloc[:, 0], ref_data.iloc[:, 9], label=f'MRI', linewidth=5, linestyle='--', color='gray')

    # 軸ラベルと設定
    plt.xlabel('$x$ [m]', fontsize=26)
    plt.ylabel('Velocity ($x$) [m/s]', fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.ylim(-0.1,0.1)
    plt.xlim(0.0, 0.05)
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


for i, (vel_file, data_file, direct_file) in enumerate(zip(data_files, data_ref_files, direct_files)):
    vel_data = pd.read_csv(vel_file)
    ref_data = pd.read_csv(data_file)
    direct_data = pd.read_csv(direct_file)

    # プロット設定
    plt.figure(figsize=(10, 6))
    
    # 各フェーズのプロット
    plt.plot(vel_data.iloc[:, 0], vel_data.iloc[:, 4], label=f'Assimilated', linewidth=4, linestyle='-', color='black')
    #plt.plot(vel_data.iloc[:, 0], direct_data.iloc[:, 4], label=f'Direct (Assimilated BC)', linewidth=5, linestyle='-', color='blue')
    plt.plot(ref_data.iloc[:, 0], ref_data.iloc[:, 10], label=f'MRI', linewidth=4, linestyle='--', color='gray')

    # 軸ラベルと設定
    plt.xlabel('$x$ [m]', fontsize=26)
    plt.ylabel('Velocity ($y$) [m/s]', fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.ylim(0.0,0.7)
    plt.xlim(0.0, 0.05)
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

for i, (vel_file, data_file, direct_file) in enumerate(zip(data_files, data_ref_files, direct_files)):
    vel_data = pd.read_csv(vel_file)
    ref_data = pd.read_csv(data_file)
    direct_data = pd.read_csv(direct_file)

    # プロット設定
    plt.figure(figsize=(10, 6))
    
    # 各フェーズのプロット
    plt.plot(vel_data.iloc[:, 0], vel_data.iloc[:, 5], label=f'Assimilated', linewidth=5, linestyle='-', color='black')
    #plt.plot(vel_data.iloc[:, 0], direct_data.iloc[:, 5], label=f'Direct (Assimilated BC)', linewidth=5, linestyle='-', color='blue')
    plt.plot(ref_data.iloc[:, 0], ref_data.iloc[:, 11], label=f'MRI', linewidth=5, linestyle='--', color='gray')

    # 軸ラベルと設定
    plt.xlabel('$x$ [m]', fontsize=26)
    plt.ylabel('Velocity ($z$) [m/s]', fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.ylim(-0.1,0.1)
    plt.xlim(0.0, 0.05)
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

