import struct
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'ieee'])

def read_bin_file(filename, data_type='d'): 
    with open(filename, 'rb') as f:
        width = struct.unpack('i', f.read(4))[0]
        height = struct.unpack('i', f.read(4))[0]
        
        num_elements = width * height
        data_size = struct.calcsize(data_type)
        
        data = struct.unpack(f'{num_elements}{data_type}', f.read(num_elements * data_size))

    return width, height, data

def calculate_mae(data_opt, data_ref):
    maes = []
    for (width_opt, height_opt, data_opt_vals), (width_ref, height_ref, data_ref_vals) in zip(data_opt, data_ref):
        data_opt_array = np.array(data_opt_vals)
        data_ref_array = np.array(data_ref_vals)
        
        mae = np.mean(np.abs(data_opt_array - data_ref_array))
        maes.append(mae)
    return maes

def process_folder(opt_folder, ref_folder, ref_pattern='velocityReference_*.bin'):
    # Optimized folder files
    file_pattern = os.path.join(opt_folder, 'velocity_*.bin')
    file_list_opt = sorted(
        [f for f in glob.glob(file_pattern)],
        key=lambda x: -1 if 'initial' in os.path.basename(x) else int(os.path.splitext(os.path.basename(x))[0].split('_')[-1])
    )

    all_opt_data = []
    for filename in file_list_opt:
        width, height, data = read_bin_file(filename)
        all_opt_data.append((width, height, data))  

    file_pattern = os.path.join(ref_folder, ref_pattern)
    file_list_ref = sorted(
        [f for f in glob.glob(file_pattern)],
        key=lambda x: -1 if 'initial' in os.path.basename(x) else int(os.path.splitext(os.path.basename(x))[0].split('_')[-1])
    )

    all_ref_data = []
    for filename in file_list_ref:
        width, height, data = read_bin_file(filename)
        all_ref_data.append((width, height, data))  
    
    return calculate_mae(all_opt_data, all_ref_data)

def plot_mae(mae_values, title, filename):
    plt.figure()
    plt.plot(mae_values, marker='o', linestyle='-', label=title)
    plt.xlabel('Time Step')
    plt.ylabel('Mean Absolute Error (MAE)')
    plt.title(f'MAE over Time - {title}')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

def plot_mae_comparison(mae_values_24x16x16, mae_values_12x8x8, mae_values_36x24x24, filename):
    plt.figure(figsize=(5, 2.5))
    x_values_36x24x24 = range(len(mae_values_36x24x24))
    x_values_24x16x16 = range(len(mae_values_24x16x16))
    x_values_12x8x8 = range(len(mae_values_12x8x8))
    
    plt.plot(x_values_36x24x24, mae_values_36x24x24, 
             marker='o', markersize=3.5, markeredgecolor='black', markeredgewidth=0.5, 
             markerfacecolor='black', linestyle='-', color='black', linewidth=0.7, label='Normalized voxel size: 0.042')
    plt.plot(x_values_24x16x16, mae_values_24x16x16, 
             marker='o', markersize=3.5, markeredgecolor='black', markeredgewidth=0.5, 
             markerfacecolor='gray', linestyle='-', color='black', linewidth=0.7, label='Normalized voxel size: 0.063')
    plt.plot(x_values_12x8x8, mae_values_12x8x8, 
             marker='o', markersize=3.5, markeredgecolor='black', markeredgewidth=0.5, 
             markerfacecolor='white', linestyle='-', color='black', linewidth=0.7, label='Normalized voxel size: 0.125')
    
    plt.xlabel('Time Step[-]')
    plt.ylabel('MAE[-]')
    plt.ylim(top=0.035, bottom=0.0)
    
    legend = plt.legend(loc='upper right', frameon=True)
    legend.get_frame().set_edgecolor('black')  
    legend.get_frame().set_facecolor('white')   
    legend.get_frame().set_linewidth(0.5)      

    plt.grid(True)
    plt.savefig(filename)
    plt.close()

opt_folder_36x24x24 = '../4dvar/output/stenosis_ds_ave_36x24x24/optimized'
ref_folder_36x24x24 = '../../direct/voxelDataCreation/output/stenosis_ds_12x8x8/bin'
mae_values_36x24x24 = process_folder(opt_folder_36x24x24, ref_folder_36x24x24)

opt_folder_24x16x16 = '../4dvar/output/stenosis_ds_ave_24x16x16/optimized'
ref_folder_24x16x16 = '../../direct/voxelDataCreation/output/stenosis_ds_12x8x8/bin'
mae_values_24x16x16 = process_folder(opt_folder_24x16x16, ref_folder_24x16x16)

opt_folder_12x8x8 = '../4dvar/output/stenosis_ds_ave_12x8x8/optimized'
ref_folder_12x8x8 = '../../direct/voxelDataCreation/output/stenosis_ds_12x8x8/bin'
mae_values_12x8x8 = process_folder(opt_folder_12x8x8, ref_folder_12x8x8)

for i, mae in enumerate(mae_values_36x24x24):
    print(f"Time Step {i} (36x24x24): MAE = {mae}")
for i, mae in enumerate(mae_values_24x16x16):
    print(f"Time Step {i} (4416x16): MAE = {mae}")
for i, mae in enumerate(mae_values_12x8x8):
    print(f"Time Step {i} (12x8x8): MAE = {mae}")

plot_mae_comparison(mae_values_24x16x16, mae_values_12x8x8, mae_values_36x24x24, 'ds_ave_mae_comparison.png')