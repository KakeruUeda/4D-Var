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
    file_pattern_opt = os.path.join(opt_folder, 'velocity_*.bin')
    file_list_opt = sorted(
        [f for f in glob.glob(file_pattern_opt)],
        key=lambda x: -1 if 'initial' in os.path.basename(x) else int(os.path.splitext(os.path.basename(x))[0].split('_')[-1])
    )

    all_opt_data = [read_bin_file(filename) for filename in file_list_opt]

    file_pattern_ref = os.path.join(ref_folder, ref_pattern)
    file_list_ref = sorted(
        [f for f in glob.glob(file_pattern_ref)],
        key=lambda x: -1 if 'initial' in os.path.basename(x) else int(os.path.splitext(os.path.basename(x))[0].split('_')[-1])
    )

    all_ref_data = [read_bin_file(filename) for filename in file_list_ref]
    
    return calculate_mae(all_opt_data, all_ref_data)

def plot_mae_vs_resolution(mae_ave_values_dict, mae_int_values_dict, filename):
    resolutions = [0.0417, 0.0625, 0.125]

    ave_avg_maes = [np.mean(mae_values) for mae_values in mae_ave_values_dict.values()]
    int_avg_maes = [np.mean(mae_values) for mae_values in mae_int_values_dict.values()]
    
    plt.figure(figsize=(5, 2.5))
    plt.plot(resolutions, ave_avg_maes, marker='o', linestyle='-', color='black', label=r'$\mathcal{T}^{ave}$')
    plt.plot(resolutions, int_avg_maes, marker='o', linestyle='-', color='blue', label=r'$\mathcal{T}^{int}$')  

    legend = plt.legend(loc='upper left', frameon=True, fontsize=10)
    legend.get_frame().set_edgecolor('black')  
    legend.get_frame().set_facecolor('white')   
    legend.get_frame().set_linewidth(0.5)    
    
    plt.xlabel('Normalized voxel size[-]')
    plt.ylabel('Time-avereged MAE[-]')
    plt.ylim(top=0.025, bottom=0.0)

    plt.grid(True)
    plt.savefig(filename)
    plt.close()

ref_folder = '../../direct/voxelDataCreation/output/stenosis_ds_36x24x24/bin'

opt_ave_folder_36x24x24 = '../4dvar/output/stenosis_ds_ave_36x24x24/optimized'
opt_int_folder_36x24x24 = '../4dvar/output/stenosis_ds_int_36x24x24/optimized'
mae_ave_values_36x24x24 = process_folder(opt_ave_folder_36x24x24, ref_folder)
mae_int_values_36x24x24 = process_folder(opt_int_folder_36x24x24, ref_folder)

opt_ave_folder_24x16x16 = '../4dvar/output/stenosis_ds_ave_24x16x16/optimized'
opt_int_folder_24x16x16 = '../4dvar/output/stenosis_ds_int_24x16x16/optimized'
mae_ave_values_24x16x16 = process_folder(opt_ave_folder_24x16x16, ref_folder)
mae_int_values_24x16x16 = process_folder(opt_int_folder_24x16x16, ref_folder)

opt_ave_folder_12x8x8 = '../4dvar/output/stenosis_ds_ave_12x8x8/optimized'
opt_int_folder_12x8x8 = '../4dvar/output/stenosis_ds_int_12x8x8/optimized'
mae_ave_values_12x8x8 = process_folder(opt_ave_folder_12x8x8, ref_folder)
mae_int_values_12x8x8 = process_folder(opt_int_folder_12x8x8, ref_folder)

mae_ave_values_dict = {
    '$66\times24\times24$': mae_ave_values_36x24x24,
    '$44\times16\times16$': mae_ave_values_24x16x16,
    '$22\times8\times8$': mae_ave_values_12x8x8
}

mae_int_values_dict = {
    '$66\times24\times24$': mae_int_values_36x24x24,
    '$44\times16\times16$': mae_int_values_24x16x16,
    '$22\times8\times8$': mae_int_values_12x8x8
}

plot_mae_vs_resolution(mae_ave_values_dict, mae_int_values_dict, 'ds_timeave_mae_comparison.png')
