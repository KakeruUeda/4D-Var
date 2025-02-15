import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scienceplots
from scipy.stats import linregress

plt.style.use(['science', 'ieee'])
# plt.style.use(['science', 'vibrant']) 

# Define paths
mri_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  
cfd_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  

mri2_path = '../../4dvar/output/Ubend_inlet_half_space_wave_time_wave_reg1e-1/optimized/'  
cfd2_path = '../../4dvar/output/Ubend_inlet_half_space_wave_time_wave_reg1e-1/optimized/'  

output_folder = "Ubend_bend1_direct_MRI_inlet"
os.makedirs(output_folder, exist_ok=True)

def sort_files_by_numeric_part(file_list):
    return sorted(file_list, key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))

# Load files
mri_files = sort_files_by_numeric_part(glob.glob(os.path.join(mri_path, 'v_mri_fluid_*.dat')))
v_mri_data = [np.loadtxt(file) for file in mri_files]

mri2_files = sort_files_by_numeric_part(glob.glob(os.path.join(mri2_path, 'v_mri_fluid_*.dat')))
v_mri2_data = [np.loadtxt(file) for file in mri2_files]

cfd_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd_path, 'v_cfd_fluid_*.dat')))
v_cfd_data = [np.loadtxt(file) for file in cfd_files]

cfd2_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd2_path, 'v_cfd_fluid_*.dat')))
v_cfd2_data = [np.loadtxt(file) for file in cfd2_files]

errors_over_time = []
errors2_over_time = []

# Process data
for time_idx in range(len(mri_files)):
    print(f"Processing Time Step: {time_idx}")

    # MRI data
    u_mri = v_mri_data[time_idx][:, 0]  
    u_mri2 = v_mri2_data[time_idx][:, 0]  

    v_mri = v_mri_data[time_idx][:, 1] 
    v_mri2 = v_mri2_data[time_idx][:, 1] 
    
    w_mri = v_mri_data[time_idx][:, 2] 
    w_mri2 = v_mri2_data[time_idx][:, 2] 

    # CFD data
    u_cfd = v_cfd_data[time_idx][:, 0] 
    u_cfd2 = v_cfd2_data[time_idx][:, 0]

    # CFD data
    v_cfd = v_cfd_data[time_idx][:, 1] 
    v_cfd2 = v_cfd2_data[time_idx][:, 1]

    # CFD data
    w_cfd = v_cfd_data[time_idx][:, 2] 
    w_cfd2 = v_cfd2_data[time_idx][:, 2]

    u_error = u_cfd - u_mri
    v_error = v_cfd - v_mri
    w_error = w_cfd - w_mri

    u_error2 = u_cfd2 - u_mri2
    v_error2 = v_cfd2 - v_mri2
    w_error2 = w_cfd2 - w_mri2

    error_mag = u_error**2 + v_error**2 + w_error**2
    error_mag2 = u_error2**2 + v_error2**2 + w_error2**2

    error_mag_sum = np.sum(error_mag)
    error_mag2_sum = np.sum(error_mag2)

    # vel_mag_mri = np.sqrt(u_mri**2 + v_mri**2 + w_mri**2)
    # vel_mag_mri_sum = np.sum(vel_mag_mri)

    # vel_mag2_mri = np.sqrt(u_mri2**2 + v_mri2**2 + w_mri2**2)
    # vel_mag2_mri_sum = np.sum(vel_mag_mri)

    # # Average errors for the current time step
    # avg_error = error_mag_sum / vel_mag_mri_sum
    # avg_error2 = error_mag2_sum / vel_mag2_mri_sum

    error_mag_sum_sqrt = np.sqrt(error_mag_sum)
    error_mag2_sum_sqrt = np.sqrt(error_mag2_sum)

    vel_mag_mri_sum = np.sum(u_mri**2 + v_mri**2 + w_mri**2)
    vel_mag_mri_sum_sqrt = np.sqrt(vel_mag_mri_sum)

    vel_mag2_mri_sum = np.sum(u_mri2**2 + v_mri2**2 + w_mri2**2)
    vel_mag2_mri_sum_sqrt = np.sqrt(vel_mag2_mri_sum)

    error = error_mag_sum_sqrt / vel_mag_mri_sum_sqrt
    error2 = error_mag2_sum_sqrt / vel_mag2_mri_sum_sqrt

    errors_over_time.append(error)
    errors2_over_time.append(error2)

dt = 0.02947812

time_values = [dt * i for i in range(len(errors_over_time))]

# Plot errors over time
plt.figure(figsize=(8, 6))
plt.plot(time_values, errors2_over_time, marker='o', linestyle = '-', linewidth = 4, color='blue', label='Inlet section')
plt.plot(time_values, errors_over_time, marker='o', linestyle = '-', linewidth = 4, color='skyblue', label='Bend section')

plt.xlabel('Time [$\mathrm{s}$]', fontsize=20)
plt.ylabel('NRMSE [-]', fontsize=20)

# 軸メモリのサイズを調整
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.legend(fontsize=20)
# plt.tight_layout()

# Save the plot
output_file = os.path.join(output_folder, "error_inlet_vs_bend.png")
plt.savefig(output_file)
print(f"Plot saved to {output_file}")