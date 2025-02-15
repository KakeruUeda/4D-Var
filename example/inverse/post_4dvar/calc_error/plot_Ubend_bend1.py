import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scienceplots
from scipy.stats import linregress

plt.style.use(['science', 'ieee'])

# Define paths
mri_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  
cfd_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  
cfd2_path = '../../4dvar/output/DataCreation_bend1_MRI_flowRate/optimized/'  
cfd3_path = '../../4dvar/output/DataCreation_bend1_assimilated/optimized/'  

# 追加: maskのパス
# mask_path = '../../4dvar/input/Ubend/bend1/data_half/'
mask_path = 'mask/'

output_folder = "Ubend_bend1_direct_MRI_inlet"
os.makedirs(output_folder, exist_ok=True)

def sort_files_by_numeric_part(file_list):
    return sorted(file_list, key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))

#--- (1) ファイル一覧の取得 ---#
mri_files = sort_files_by_numeric_part(glob.glob(os.path.join(mri_path, 'v_mri_fluid_*.dat')))
v_mri_data = [np.loadtxt(file) for file in mri_files]

cfd_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd_path, 'v_cfd_fluid_*.dat')))
v_cfd_data = [np.loadtxt(file) for file in cfd_files]

cfd2_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd2_path, 'v_cfd_fluid_*.dat')))
v_cfd2_data = [np.loadtxt(file) for file in cfd2_files]

cfd3_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd3_path, 'v_cfd_fluid_*.dat')))
v_cfd3_data = [np.loadtxt(file) for file in cfd3_files]

#--- (2) mask の読み込み ---#
# 例えば "mask.dat" というファイル名なら以下のように読み込む
mask_file = os.path.join(mask_path, 'mask_x.dat')
mask = np.loadtxt(mask_file)

# mask < 0 となる要素の真偽値を持つ配列を作成
valid_indices = (mask > 0.9999)

errors_over_time  = []
errors2_over_time = []
errors3_over_time = []

#--- (3) データ処理 ---#
for time_idx in range(len(mri_files)):

    #--- (3.1) MRI data 読み込み ---#
    u_mri = v_mri_data[time_idx][:, 0]  
    v_mri = v_mri_data[time_idx][:, 1] 
    w_mri = v_mri_data[time_idx][:, 2] 

    #--- (3.2) CFD data 読み込み ---#
    u_cfd  = v_cfd_data[time_idx][:, 0] 
    v_cfd  = v_cfd_data[time_idx][:, 1] 
    w_cfd  = v_cfd_data[time_idx][:, 2] 

    u_cfd2 = v_cfd2_data[time_idx][:, 0]
    v_cfd2 = v_cfd2_data[time_idx][:, 1]
    w_cfd2 = v_cfd2_data[time_idx][:, 2]

    u_cfd3 = v_cfd3_data[time_idx][:, 0]
    v_cfd3 = v_cfd3_data[time_idx][:, 1]
    w_cfd3 = v_cfd3_data[time_idx][:, 2]

    #--- (3.3) mask < 0 の要素のみ抽出 ---#
    u_mri_masked  = u_mri[valid_indices]
    v_mri_masked  = v_mri[valid_indices]
    w_mri_masked  = w_mri[valid_indices]

    u_cfd_masked  = u_cfd[valid_indices]
    v_cfd_masked  = v_cfd[valid_indices]
    w_cfd_masked  = w_cfd[valid_indices]

    u_cfd2_masked = u_cfd2[valid_indices]
    v_cfd2_masked = v_cfd2[valid_indices]
    w_cfd2_masked = w_cfd2[valid_indices]

    u_cfd3_masked = u_cfd3[valid_indices]
    v_cfd3_masked = v_cfd3[valid_indices]
    w_cfd3_masked = w_cfd3[valid_indices]

    #--- (3.4) 誤差計算 ---#
    # (DA)
    u_error  = u_cfd_masked  - u_mri_masked
    v_error  = v_cfd_masked  - v_mri_masked
    w_error  = w_cfd_masked  - w_mri_masked
    error_mag  = u_error**2 + v_error**2 + w_error**2
    error_mag_sum = np.sum(error_mag)
    error_mag_sum_sqrt = np.sqrt(error_mag_sum)

    # (Direct (Poiseuille BC) - 参考例)
    u_error2 = u_cfd2_masked - u_mri_masked
    v_error2 = v_cfd2_masked - v_mri_masked
    w_error2 = w_cfd2_masked - w_mri_masked
    error_mag2 = u_error2**2 + v_error2**2 + w_error2**2
    error_mag2_sum = np.sum(error_mag2)
    error_mag2_sum_sqrt = np.sqrt(error_mag2_sum)

    # (Direct (Assimilated BC))
    u_error3 = u_cfd3_masked - u_mri_masked
    v_error3 = v_cfd3_masked - v_mri_masked
    w_error3 = w_cfd3_masked - w_mri_masked
    error_mag3 = u_error3**2 + v_error3**2 + w_error3**2
    error_mag3_sum = np.sum(error_mag3)
    error_mag3_sum_sqrt = np.sqrt(error_mag3_sum)

    #--- (3.5) MRI速度ベクトルの大きさの和 (NRMSE用) ---#
    vel_mag_mri = u_mri_masked**2 + v_mri_masked**2 + w_mri_masked**2
    vel_mag_mri_sum = np.sum(vel_mag_mri)
    vel_mag_mri_sum_sqrt = np.sqrt(vel_mag_mri_sum)

    #--- (3.6) 正規化誤差(NRMSE) ---#
    avg_error  = error_mag_sum_sqrt  / vel_mag_mri_sum_sqrt
    avg_error2 = error_mag2_sum_sqrt / vel_mag_mri_sum_sqrt
    avg_error3 = error_mag3_sum_sqrt / vel_mag_mri_sum_sqrt

    errors_over_time.append(avg_error)
    errors2_over_time.append(avg_error2)
    errors3_over_time.append(avg_error3)

#--- (4) 時系列プロット ---#
dt = 0.02947812
time_values = [dt * i for i in range(len(errors_over_time))]

plt.figure(figsize=(8, 6))
plt.plot(time_values, errors_over_time,  marker='o', linestyle='-', color='black', linewidth=4, label='DA')
plt.plot(time_values, errors3_over_time, marker='o', linestyle='-', color='blue',  linewidth=4, label='Direct')
# 必要に応じて別の BC の曲線を再度描画
# plt.plot(time_values, errors2_over_time, marker='o', linestyle='-', color='lightgreen', linewidth=4, label='Direct ($\\mathbf{V}^{\\mathrm{p}}$)')

plt.xlabel('Time [$\\mathrm{s}$]', fontsize=20)
plt.ylabel('NRMSE [-]', fontsize=20)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.legend(fontsize=20)

output_file = os.path.join(output_folder, "error_over_time.png")
plt.savefig(output_file)
print(f"Plot saved to {output_file}")
