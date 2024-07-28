import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'ieee'])

file_path_ref = r'output/flowRate8x8x8_interpolation/flowRateVelRef.dat'
file_path_opt_ave = r'output/flowRate8x8x8_average/flowRateVelOpt.dat'
file_path_opt_int = r'output/flowRate8x8x8_interpolation/flowRateVelOpt.dat'

save_path = r'flowRate8x8x8.png'

data_ref = np.loadtxt(file_path_ref)
data_opt_ave = np.loadtxt(file_path_opt_ave)
data_opt_int = np.loadtxt(file_path_opt_int)

print(data_ref[:5])
print(data_opt_ave[:5])
print(data_opt_int[:5])

T = max(data_ref[:, 0].max(), data_opt_ave[:, 0].max())

plt.figure(figsize=(7, 5))
plt.scatter(data_ref[:, 0], data_ref[:, 1], label='Reference', marker='o', s=70)
plt.scatter(data_opt_ave[:, 0], data_opt_ave[:, 1], label='Smoothing', marker='x', s=70, linewidths=2)
plt.scatter(data_opt_int[:, 0], data_opt_int[:, 1], label='Interpolation', marker='x', s=70, linewidths=2)

plt.xlabel(r'Time [-]', fontsize=18)
plt.ylabel(r'Flow rate [-]', fontsize=18)

ticks = np.linspace(0, T, 11)
tick_labels = [f'{i/10:.1f}T' for i in range(11)]
plt.xticks(ticks, tick_labels, fontsize=14)

plt.yticks(fontsize=18)

plt.legend(fontsize=18)

plt.grid(True)

plt.savefig(save_path, dpi=300)
