#include "../../../../src/lib/common/Array.h"
#include <array>
#include <cmath>
#include <filesystem>  // C++17以上で使用可能
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>

namespace fs = std::filesystem;

int main()
{
  // ディレクトリ名やファイル名
  std::string results_dir = "../../../inverse/4dvar/output/Ubend_inlet_space_wave_time_wave_reg1e-1/optimized/";

  //==============================
  // 3. MRI 側流速データ読み込み
  //==============================
  std::vector<std::vector<std::vector<double>>> data;  
  // MRI 側のファイルは "v_mri_fluid_0.dat", "v_mri_fluid_1.dat", ...
  for(int t = 0;; t++) {
    std::string file_name = results_dir + "v_mri_fluid_" + std::to_string(t) + ".dat";

    // ファイルが存在しないときはループ終了
    if(!fs::exists(file_name)) {
      break;
    }
    std::ifstream file(file_name);
    if(!file.is_open()) {
      std::cerr << "Failed to open file: " << file_name << std::endl;
      continue;
    }

    std::vector<std::vector<double>> plane;
    std::string line;
    while(std::getline(file, line)) {
      std::stringstream ss(line);
      std::vector<double> velocities;
      double value;
      while(ss >> value) {
        velocities.push_back(value);
      }
      plane.push_back(velocities);
    }
    data.push_back(plane);
  }

  //==============================
  // 4. 格子情報などの設定
  //==============================
  // MRI データの格子情報
  std::array<int, 3> nxData = {27, 27, 27};
  std::array<double, 3> lxData = {0.050625, 0.050625, 0.0513};
  std::array<double, 3> dxData = {
    lxData[0] / nxData[0],
    lxData[1] / nxData[1],
    lxData[2] / nxData[2]
  };

  double dt_mri = 0.02947812;

  // 断面位置(例として中心 j=0 としているが、本来は適宜指定)
  int j_center_data = 0;
  int j_center_cfd = 0;

  //==============================
  // 5. MRI 側の流量(Qmri)を計算
  //==============================
  std::vector<double> flow_rate_mri(data.size(), 0.0);

  for(size_t t = 0; t < data.size(); ++t) {
    double flow_rate_data = 0.0;

    for(int k = 0; k < nxData[2]; k++) {
      for(int j = 0; j < nxData[1]; j++) {
        for(int i = 0; i < nxData[0]; i++) {
          if(j == j_center_data) {
            int index = k * nxData[0] * nxData[1] + j * nxData[0] + i;
            flow_rate_data += data[t][index][1] * dxData[0] * dxData[2];
          }
        }
      }
    }
    flow_rate_mri[t] = flow_rate_data;
  }

  //==============================
  // 7. 結果をファイル出力(流量そのもの)
  //==============================
  std::ofstream flow_rate_mri_file("flow_rate_mri.dat");
  
  // MRI側(時刻 t*dt_mri)の流量
  for(size_t t = 0; t < flow_rate_mri.size(); ++t) {
    double time_mri = t * dt_mri;
    flow_rate_mri_file << time_mri << " " << flow_rate_mri[t] << "\n";
  }
  flow_rate_mri_file.close();

  //==============================
  // 9. 正常終了
  //==============================
  return 0;
}