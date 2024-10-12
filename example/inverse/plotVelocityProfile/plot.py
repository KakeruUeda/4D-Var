import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'ieee'])

#vel1 = pd.read_csv('../../../../..//Results/4D-Var/Inverse/output/Ubend/average/velocity_noreg/velocity_profile/velocity_xy_y=0_step0.csv')
vel2 = pd.read_csv('../../../../..//Results/4D-Var/Inverse/output/Ubend/average/velocity_profile/velocity_xy_y=0.0253125_step0.csv')
data = pd.read_csv('../../../../..//Results/4D-Var/Inverse/output/Ubend/average/velocity_profile/data_xy_y=0.0253125_phase0.csv')

save_path = r'velocityProfile_y=0.0253125_step0.png'

plt.figure(figsize=(10, 6))
#plt.plot(vel1.iloc[:, 0], vel1.iloc[:, 4], label=r'Optimized', linewidth=2.7, linestyle='-', color='red')
plt.plot(vel2.iloc[:, 0], vel2.iloc[:, 4], label=r'Optimized', linewidth=2.7, linestyle='-', color='blue')
plt.plot(data.iloc[:, 0], data.iloc[:, 4], label='Data', linewidth=2.7, linestyle='--', color='darkblue')

plt.xlabel('x coordinates [m]', fontsize=20)
plt.ylabel('velocity magnitude [m/s]', fontsize=20)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.ylim(top=0.5)

plt.legend(fontsize=20, loc='upper right', frameon=True)

plt.grid(True)

plt.savefig(save_path, dpi=300)
