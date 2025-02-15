import matplotlib.pyplot as plt
import numpy as np
import scienceplots
from matplotlib.ticker import ScalarFormatter

plt.style.use(['science', 'ieee'])

# ファイルからデータを読み込む関数
def read_flow_rate(file_path):
    time = []
    flow_rate = []
    with open(file_path, 'r') as file:
        for line in file:
            t, rate = map(float, line.strip().split())
            time.append(t)
            flow_rate.append(rate)
    return time, flow_rate

# ファイルパス
flow_rate_mri_file = "flow_rate_bend_mri.dat"
flow_rate_direct_assimilated_BC_file = "flow_rate_bend_direct_assimilated_BC.dat"
flow_rate_direct_MRI_BC_file = "flow_rate_bend_direct_MRI_BC.dat"
flow_rate_DA_file = "flow_rate_bend_assimilated.dat"

# データの読み込み
time_mri, flow_rate_mri = read_flow_rate(flow_rate_mri_file)
time_direct, flow_rate_direct_assimilated_BC = read_flow_rate(flow_rate_direct_assimilated_BC_file)
time_direct, flow_rate_direct_MRI_BC = read_flow_rate(flow_rate_direct_MRI_BC_file)
time_direct, flow_rate_DA = read_flow_rate(flow_rate_DA_file)

# 流量を mm 単位に変換
# flow_rate_mri = np.array(flow_rate_mri) 
# flow_rate_direct_assimilated_BC = np.array(flow_rate_direct_assimilated_BC) 
# flow_rate_direct_MRI_BC = np.array(flow_rate_direct_MRI_BC) 

plt.figure(figsize=(12, 8))

# MRIプロットに丸マーカーを追加
plt.plot(time_mri, flow_rate_mri, label="MRI", linewidth=5, color='gray', linestyle='--', marker='o', markersize=10)
plt.plot(time_direct, flow_rate_DA, label='DA', linewidth=5, linestyle='-', color='black')
plt.plot(time_direct, flow_rate_direct_assimilated_BC, label='Direct ($\mathbf{V}^{\mathrm{a}}$)', linewidth=5, linestyle='-', color='blue')
# plt.plot(time_direct, flow_rate_direct_MRI_BC, label='Direct ($\mathbf{V}^{\mathrm{p}}$)', linewidth=5, linestyle='-', color='lightgreen')

# 軸の設定
ax = plt.gca()

# # 外枠（top, right）を非表示
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# # x軸とy軸を太くする
# ax.spines['bottom'].set_linewidth(2.5)
# ax.spines['left'].set_linewidth(2.5)

# 軸のラベルとタイトル
plt.xlabel("Time $[\mathrm{s}]$", fontsize=28)
plt.ylabel("Flow rate $[\mathrm{m}^3/\mathrm{s}]$", fontsize=28)

# 軸メモリのサイズを調整
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)

# Y軸範囲設定
plt.ylim(0, 0.0002)

# Y軸の指数表記設定
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
ax.yaxis.offsetText.set_fontsize(28)  # ここで指数部分のフォントサイズを変更

# 凡例
plt.legend(fontsize=28, loc='upper right')

# 保存
save_path = r'flow_rate_bend_direct_assimilated_BC.png'
plt.savefig(save_path, dpi=300)
