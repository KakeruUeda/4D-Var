import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'ieee'])

file_path = '../4dvar/output/Ubend_ave_reg2e-4/dat/costFunction.dat'
data = np.loadtxt(file_path)

save_path = r'costFunction_average.png'

plt.figure(figsize=(10, 6))

plt.plot(range(data.shape[0]), data[:, 7], label='Iteration history of cost function', linewidth=2.7, linestyle='-')
#plt.plot(range(data.shape[0]), data[:, 2], label='term 3', linewidth=2.7, linestyle='-')

plt.xlabel(r'Iteration [-]', fontsize=24)
plt.ylabel(r'Value [-]', fontsize=24)

plt.yscale('log')

plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.legend(fontsize=20)
plt.grid(True)

plt.savefig(save_path, dpi=300)
