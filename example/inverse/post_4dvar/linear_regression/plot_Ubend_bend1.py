import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scienceplots
from scipy.stats import linregress

plt.style.use(['science', 'ieee'])

mri_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  
cfd_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  

# mri_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  
# cfd_path = '../../4dvar/output/DataCreation_bend1_assimilated/optimized/'  

# mri_path = '../../4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/'  
# cfd_path = '../../4dvar/output/DataCreation_bend1_MRI_flowRate/optimized/'  

output_folder = "Ubend_bend1_DA"
os.makedirs(output_folder, exist_ok=True)

def sort_files_by_numeric_part(file_list):
    return sorted(file_list, key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))

mri_files = sort_files_by_numeric_part(glob.glob(os.path.join(mri_path, 'v_mri_fluid_*.dat')))
v_mri_data = [np.loadtxt(file) for file in mri_files]

cfd_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd_path, 'v_cfd_fluid_*.dat')))
v_cfd_data = [np.loadtxt(file) for file in cfd_files]

regression_coefficients = []
r_squared_values = []

for time_idx, (mri_file, cfd_file) in enumerate(zip(mri_files, cfd_files)):
    print(f"Processing Time Step: {time_idx}")
    print(f"  MRI Data File: {mri_file}")
    print(f"  CFD Data File: {cfd_file}")

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
    # --- プロット ---
    plt.figure(figsize=(12, 10))
    plt.scatter(velocity_mri, velocity_cfd, marker='+', alpha=0.5, s=200, color='red')

    plt.xlabel(r'MRI velocity ($x$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=44)
    plt.ylabel(r'Assimilated velocity ($x$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=44)

    plt.xlim(-0.3, 0.7)
    plt.ylim(-0.3, 0.7)
    
    plt.xticks(fontsize=44, fontweight='bold')
    plt.yticks(fontsize=44, fontweight='bold')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    ax.tick_params(
        direction='in',      
        length=10,         
        width=3,         
        pad=8            
    )

    # --- 回帰分析 ---
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
        fontsize=44, 
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

    ax = plt.gca()
    ax.set_aspect(1.0) 

    ax.xaxis.set_major_locator(plt.MaxNLocator(5))  
    ax.yaxis.set_major_locator(plt.MaxNLocator(5)) 

    plt.tight_layout()

    plt.legend(fontsize=44, loc='upper left')

    output_path = os.path.join(output_folder, f"u_linear_regression_time{time_idx}.png")
    plt.savefig(output_path)
    plt.close()
    print(f"  Plot saved to {output_path}")

# --- 平均スロープとR^2を計算 ---
average_slope = np.mean(regression_coefficients)
average_r_squared = np.mean(r_squared_values)

print(f"Average Regression Coefficient (Slope): {average_slope:.3f}")
print(f"Average R^2 Value: {average_r_squared:.3f}")



# ---- y dir. ------
for time_idx, (mri_file, cfd_file) in enumerate(zip(mri_files, cfd_files)):
    print(f"Processing Time Step: {time_idx}")
    print(f"  MRI Data File: {mri_file}")
    print(f"  CFD Data File: {cfd_file}")

    u_mri = v_mri_data[time_idx][:, 0]  
    v_mri = v_mri_data[time_idx][:, 1] 
    w_mri = v_mri_data[time_idx][:, 2] 

    u_cfd = v_cfd_data[time_idx][:, 0] 
    v_cfd = v_cfd_data[time_idx][:, 1] 
    w_cfd = v_cfd_data[time_idx][:, 2] 

    v_mag_mri = np.sqrt(u_mri**2 + v_mri**2 + w_mri**2) 
    v_mag_cfd = np.sqrt(u_cfd**2 + v_cfd**2 + w_cfd**2)

    velocity_mri = v_mri
    velocity_cfd = v_cfd
    # --- プロット ---
    plt.figure(figsize=(12, 10))
    plt.scatter(velocity_mri, velocity_cfd, marker='+', alpha=0.5, s=200, color='blue')

    plt.xlabel(r'MRI velocity ($y$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=44)
    plt.ylabel(r'Assimilated velocity ($y$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=44)

    plt.xlim(-0.3, 0.7)
    plt.ylim(-0.3, 0.7)
    
    plt.xticks(fontsize=44, fontweight='bold')
    plt.yticks(fontsize=44, fontweight='bold')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    ax.tick_params(
        direction='in',      
        length=10,         
        width=3,         
        pad=8            
    )

    # --- 回帰分析 ---
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
        fontsize=44, 
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

    ax = plt.gca()
    ax.set_aspect(1.0) 

    ax.xaxis.set_major_locator(plt.MaxNLocator(5))  
    ax.yaxis.set_major_locator(plt.MaxNLocator(5)) 

    plt.tight_layout()

    plt.legend(fontsize=44, loc='upper left')

    output_path = os.path.join(output_folder, f"v_linear_regression_time{time_idx}.png")
    plt.savefig(output_path)
    plt.close()
    print(f"  Plot saved to {output_path}")


# ---- y dir. ------
for time_idx, (mri_file, cfd_file) in enumerate(zip(mri_files, cfd_files)):
    print(f"Processing Time Step: {time_idx}")
    print(f"  MRI Data File: {mri_file}")
    print(f"  CFD Data File: {cfd_file}")

    u_mri = v_mri_data[time_idx][:, 0]  
    v_mri = v_mri_data[time_idx][:, 1] 
    w_mri = v_mri_data[time_idx][:, 2] 

    u_cfd = v_cfd_data[time_idx][:, 0] 
    v_cfd = v_cfd_data[time_idx][:, 1] 
    w_cfd = v_cfd_data[time_idx][:, 2] 

    v_mag_mri = np.sqrt(u_mri**2 + v_mri**2 + w_mri**2) 
    v_mag_cfd = np.sqrt(u_cfd**2 + v_cfd**2 + w_cfd**2)

    velocity_mri = w_mri
    velocity_cfd = w_cfd

    # --- プロット ---
    plt.figure(figsize=(12, 10))
    plt.scatter(velocity_mri, velocity_cfd, marker='+', alpha=0.5, s=200, color='lightgreen')

    plt.xlabel(r'MRI velocity ($z$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=44)
    plt.ylabel(r'Assimilated velocity ($z$) $\ [\mathrm{m}/\mathrm{s}]$', fontsize=44)

    plt.xlim(-0.3, 0.3)
    plt.ylim(-0.3, 0.3)
    
    plt.xticks(fontsize=44, fontweight='bold')
    plt.yticks(fontsize=44, fontweight='bold')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    ax.tick_params(
        direction='in',      
        length=10,         
        width=3,         
        pad=8            
    )

    # --- 回帰分析 ---
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
        fontsize=44, 
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

    ax = plt.gca()
    ax.set_aspect(1.0) 

    ax.xaxis.set_major_locator(plt.MaxNLocator(5))  
    ax.yaxis.set_major_locator(plt.MaxNLocator(5)) 

    plt.tight_layout()

    plt.legend(fontsize=36, loc='upper left')

    output_path = os.path.join(output_folder, f"w_linear_regression_time{time_idx}.png")
    plt.savefig(output_path)
    plt.close()
    print(f"  Plot saved to {output_path}")