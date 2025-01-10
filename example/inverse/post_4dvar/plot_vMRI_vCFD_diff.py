import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scienceplots

plt.style.use(['science', 'ieee'])

mri_path = '../4dvar/output/Ubend_ave_reg2e-4/optimized/'
cfd_path = '../4dvar/output/Ubend_ave_reg2e-4/optimized/' 

output_folder = "regression_v_diff"
os.makedirs(output_folder, exist_ok=True)

def sort_files_by_numeric_part(file_list):
    return sorted(file_list, key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split('_')[-1]))

mri_files = sort_files_by_numeric_part(glob.glob(os.path.join(mri_path, 'v_mri_*.dat')))
v_mri_data = [np.loadtxt(file) for file in mri_files]  

cfd_files = sort_files_by_numeric_part(glob.glob(os.path.join(cfd_path, 'v_cfd_*.dat')))
v_cfd_data = [np.loadtxt(file) for file in cfd_files]  

for time_idx, (mri_file, cfd_file) in enumerate(zip(mri_files, cfd_files)):
    print(f"Processing Time Step: {time_idx}")
    print(f"  MRI Data File: {mri_file}")
    print(f"  CFD Data File: {cfd_file}")

    column_index = 1

    v_mri = v_mri_data[time_idx][:, column_index] 
    v_cfd = v_cfd_data[time_idx][:, column_index] 

    v_diff = v_cfd - v_mri 

    mean_diff = np.mean(v_diff)
    std_diff = np.std(v_diff)

    two_std_diff = 2 * std_diff
    minus_two_std_diff = -two_std_diff

    upper_2sd = mean_diff + 2 * std_diff
    lower_2sd = mean_diff - 2 * std_diff

    plt.figure(figsize=(12, 10))
    plt.scatter(v_mri, v_diff, marker='+', alpha=0.5, s=200, color='blue')
    plt.axhline(mean_diff, color='gray', linestyle='--', linewidth=3)
    plt.axhline(upper_2sd, color='gray', linestyle='--', linewidth=2)
    plt.axhline(lower_2sd, color='gray', linestyle='--', linewidth=2)

    plt.xlabel(r'$\mathbf{v}_\mathrm{m} \ [\mathrm{m}/\mathrm{s}]$', fontsize=40)
    plt.ylabel(r'$\mathcal{T}^{ave} \mathbf{v} - \mathbf{v}_\mathrm{m} \ [\mathrm{m}/\mathrm{s}]$', fontsize=40)

    plt.axis([-0.5, 0.5, -0.5, 0.5])

    plt.xticks(fontsize=40, fontweight='bold')
    plt.yticks(fontsize=40, fontweight='bold')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    plt.text(
        0.05, 0.85, 
        f"Mean±2SD: \n {mean_diff:.2f}±{two_std_diff:.2f} [m/s]",
        fontsize=40, 
        color='black', 
        transform=plt.gca().transAxes,
    )

    plt.gca().set_aspect('equal', adjustable='box')

    output_path = os.path.join(output_folder, f"linear_regression_diff_time{time_idx}.png")
    plt.savefig(output_path)
    plt.close()
    print(f"  Plot saved to {output_path}")
