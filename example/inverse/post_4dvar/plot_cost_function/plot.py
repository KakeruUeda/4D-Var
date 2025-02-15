import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'vibrant']) 

file_path = '../../4dvar/output/Ubend_inlet_half_space_wave_time_wave_reg1e-1/dat/costFunction.dat'
data = np.loadtxt(file_path)

save_path = r'costFunction_inlet.png'

plt.figure(figsize=(10, 6))

plt.plot(range(data.shape[0]), data[:, 0], label='term 1', linewidth=2.7, linestyle='-')
plt.plot(range(data.shape[0]), data[:, 3], label='term 2', linewidth=2.7, linestyle='-')
plt.plot(range(data.shape[0]), data[:, 2], label='term 3', linewidth=2.7, linestyle='-')
plt.plot(range(data.shape[0]), data[:, 6], label='term 4', linewidth=2.7, linestyle='-')

plt.xlabel(r'Iteration [-]', fontsize=24)
plt.ylabel(r'Value [-]', fontsize=24)

plt.yscale('log')

plt.ylim(top=10)

plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.legend(fontsize=20, loc='upper right', ncol=2)  
plt.grid(True)

plt.savefig(save_path, dpi=300)
