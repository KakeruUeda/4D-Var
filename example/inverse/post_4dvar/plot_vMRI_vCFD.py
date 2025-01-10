import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scienceplots
from scipy.stats import linregress

plt.style.use(['science', 'ieee'])

mri_path = '../4dvar/output/Ubend_full_ave_1e-4/optimized/'  
cfd_path = '../4dvar/output/Ubend_full_ave_1e-4/optimized/'  

output_folder = "Ubend_full_ave_regression_u"
os.makedirs(output_folder, exist_ok=True)

def sort_files_by_numeric_part(file_list):
    return sorted(file_list, key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))

mri_files = sort_files_by_numeric_part(glob.glob(os.path.join(mri_path, 'v_mri_*.dat')))
v_mri_data = [np.loadtxt(file) for file in mri_files]  

cfd_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd_path, 'v_cfd_*.dat')))
v_cfd_data = [np.loadtxt(file) for file in cfd_files] 

regression_coefficients = []
r_squared_values = []

for time_idx, (mri_file, cfd_file) in enumerate(zip(mri_files, cfd_files)):
    print(f"Processing Time Step: {time_idx}")
    print(f"  MRI Data File: {mri_file}")
    print(f"  CFD Data File: {cfd_file}")

    column_index = 0

    u_mri = v_mri_data[time_idx][:, 0]  
    v_mri = v_mri_data[time_idx][:, 1] 
    w_mri = v_mri_data[time_idx][:, 2] 

    u_cfd = v_cfd_data[time_idx][:, 0] 
    v_cfd = v_cfd_data[time_idx][:, 1] 
    w_cfd = v_cfd_data[time_idx][:, 2] 

    v_mag_mri = np.sqrt(u_mri**2 + v_mri**2 + w_mri**2) 
    v_mag_cfd = np.sqrt(u_cfd**2 + v_cfd**2 + w_cfd**2)

    velocity_mri = u_mri
    velocity_cfd = u_cfd

    plt.figure(figsize=(12, 10))
    plt.scatter(velocity_mri, velocity_cfd, marker='+', alpha=0.5, s=200, color='blue') 

    plt.xlabel(r'MRI velocity ($x$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=40)
    plt.ylabel(r'Numerical velocity ($x$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=40)

    plt.xlim(-0.25, 0.25)
    plt.ylim(-0.25, 0.25)
    
    # plt.xticks(np.arange(0, 0.6, 0.1), fontsize=40, fontweight='bold')
    # plt.yticks(np.arange(0, 0.6, 0.1), fontsize=40, fontweight='bold')

    plt.xticks(fontsize=40, fontweight='bold')
    plt.yticks(fontsize=40, fontweight='bold')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    ax.tick_params(
        direction='in',      
        length=10,         
        width=3,         
        pad=8            
    )

    slope, intercept = np.polyfit(velocity_mri, velocity_cfd, 1) 
    regression_line = slope * velocity_mri + intercept

    velocity_mri_min = np.min(velocity_mri)
    velocity_mri_max = np.max(velocity_mri)
    x_extended = np.linspace(velocity_mri_min - 100, velocity_mri_max + 100, 100)
    extended_regression_line = slope * x_extended + intercept

    r = np.corrcoef(velocity_mri, velocity_cfd)[0, 1]
    r_squared = r**2

    regression_coefficients.append(slope)
    r_squared_values.append(r_squared)

    plt.text(
        0.07, 0.80, 
        f"$r^2 = {r_squared:.2f}$", 
        fontsize=36, 
        color='black', 
        fontweight='bold',
        transform=plt.gca().transAxes 
    )

    print(f"  Time {time_idx}: r: {r:.3f}")
    print(f"  Time {time_idx}: r^2: {r_squared:.3f}")

    label = f'$y={slope:.2f}x{"+" if intercept >= 0 else ""}{intercept:.2f}$'

    plt.plot(
        x_extended,
        extended_regression_line,
        color='black',
        linestyle='--',
        linewidth=4,
        label=label
    )

    plt.tight_layout()

    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(fontsize=36, loc='upper left')

    output_path = os.path.join(output_folder, f"linear_regression_time{time_idx}.png")
    plt.savefig(output_path)
    plt.close()
    print(f"  Plot saved to {output_path}")

average_slope = np.mean(regression_coefficients)
average_r_squared = np.mean(r_squared_values)

print(f"Average Regression Coefficient (Slope): {average_slope:.3f}")
print(f"Average R^2 Value: {average_r_squared:.3f}")