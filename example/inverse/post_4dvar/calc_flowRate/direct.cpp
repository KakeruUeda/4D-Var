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
  //==============================
  // 1. パラメータや領域サイズの設定
  //==============================
  // CFDの計算結果(optimized velocity)が入るコンテナ
  std::vector<Array2D<double>> velocity_opt;
  std::vector<Array2D<double>> velocity_opt2;

  velocity_opt.resize(76);
  velocity_opt2.resize(76);

  std::string data_dir = 
    "../../../inverse/4dvar/output/Ubend_bend1_half_space_wave_time_wave_reg1e-1/optimized/";

  std::string opt_dir = 
    "../../../direct/voxelDataCreation/output/Ubend_assimilated_BC/bin/";

    std::string opt2_dir = 
    "../../../direct/voxelDataCreation/output/Ubend_MRI_BC/bin/";
    
  // MRI のファイルが何ステップあるかは実際のファイル数に依存
  // (後でファイル存在チェックで確定させる)

  // MRI の標準偏差(不要なら読み飛ばしてもOK)
  std::string std_dev_file = "std_dev.dat";
  std::vector<double> std_dev;
  {
    std::ifstream stream(std_dev_file);
    if(!stream.is_open()) {
      std::cerr << "Failed to open file: " << std_dev_file << std::endl;
      return -1;
    }
    double value;
    while(stream >> value) {
      std_dev.push_back(value);
    }
  }

  //==============================
  // 2. CFD 側流速データ読み込み
  //==============================
  for(int t = 0; t < static_cast<int>(velocity_opt.size()); t++) {
    std::string file_name = opt_dir + "velocityReference_" + std::to_string(t) + ".bin";
    velocity_opt[t].importBIN(file_name);
  }

  for(int t = 0; t < static_cast<int>(velocity_opt2.size()); t++) {
    std::string file_name = opt2_dir + "velocityReference_" + std::to_string(t) + ".bin";
    velocity_opt2[t].importBIN(file_name);
  }

  //==============================
  // 3. MRI 側流速データ読み込み
  //==============================
  std::vector<std::vector<std::vector<double>>> data;  
  // MRI 側のファイルは "v_mri_fluid_0.dat", "v_mri_fluid_1.dat", ...
  for(int t = 0;; t++) {
    std::string file_name = data_dir + "v_mri_fluid_" + std::to_string(t) + ".dat";

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
  std::array<int, 3> nxData = {40, 40, 27};
  std::array<double, 3> lxData = {0.075, 0.075, 0.0513};
  std::array<double, 3> dxData = {
    lxData[0] / nxData[0],
    lxData[1] / nxData[1],
    lxData[2] / nxData[2]
  };

  // CFD データの格子情報
  std::array<int, 3> nx = {94, 94, 64};
  std::array<double, 3> lx = {0.075, 0.075, 0.0513};
  std::array<double, 3> dx = {
    lx[0] / nx[0],
    lx[1] / nx[1],
    lx[2] / nx[2]
  };

  // 時間刻み
  double dt_cfd = 0.005895624;
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

    // data[t] は「全空間(x,y,z)」ベクトル(ただし1次元配列的に並んでいる)  
    // [k * nxData[0]*nxData[1] + j*nxData[0] + i] 行に、 velocity(例えば x,y,z) が格納されている想定
    // velocityのうち流量に使う成分: ここでは第2成分(インデックス1) = v
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
  // 6. CFD 側の流量(Qcfd)を計算
  //==============================
  std::vector<double> flow_rate_cfd(velocity_opt.size(), 0.0);

  for(size_t t = 0; t < velocity_opt.size(); ++t) {
    double flow_rate_opt = 0.0;

    for(int k = 0; k < nx[2]+1; k++) {
      for(int j = 0; j < nx[1]+1; j++) {
        for(int i = 0; i < nx[0]+1; i++) {
          if(j == j_center_cfd) {
            int index = k * (nx[0]+1) * (nx[1]+1) + j * (nx[0]+1) + i;
            flow_rate_opt += velocity_opt[t](index, 1) * dx[0] * dx[2];
          }
        }
      }
    }
    flow_rate_cfd[t] = flow_rate_opt;
  }

  std::vector<double> flow_rate_cfd2(velocity_opt2.size(), 0.0);

  for(size_t t = 0; t < velocity_opt2.size(); ++t) {
    double flow_rate_opt2 = 0.0;

    for(int k = 0; k < nx[2]+1; k++) {
      for(int j = 0; j < nx[1]+1; j++) {
        for(int i = 0; i < nx[0]+1; i++) {
          if(j == j_center_cfd) {
            int index = k * (nx[0]+1) * (nx[1]+1) + j * (nx[0]+1) + i;
            flow_rate_opt2 += velocity_opt2[t](index, 1) * dx[0] * dx[2];
          }
        }
      }
    }
    flow_rate_cfd2[t] = flow_rate_opt2;
  }

  //==============================
  // 7. 結果をファイル出力(流量そのもの)
  //==============================
  std::ofstream flow_rate_mri_file("flow_rate_bend_mri.dat");
  std::ofstream flow_rate_direct_assimilated_BC("flow_rate_bend_direct_assimilated_BC.dat");
  std::ofstream flow_rate_direct_MRI_BC("flow_rate_bend_direct_MRI_BC.dat");

  if(!flow_rate_mri_file.is_open() || !flow_rate_direct_assimilated_BC.is_open() || !flow_rate_direct_MRI_BC.is_open()) {
    std::cerr << "Failed to open output file for flow rates." << std::endl;
    return -1;
  }

  // MRI側(時刻 t*dt_mri)の流量
  for(size_t t = 0; t < flow_rate_mri.size(); ++t) {
    double time_mri = t * dt_mri;
    flow_rate_mri_file << time_mri << " " << flow_rate_mri[t] << "\n";
  }
  flow_rate_mri_file.close();

  // CFD側(時刻 t*dt_cfd)の流量
  for(size_t t = 0; t < flow_rate_cfd.size(); ++t) {
    double time_cfd = t * dt_cfd;
    flow_rate_direct_assimilated_BC << time_cfd << " " << flow_rate_cfd[t] << "\n";
  }
  flow_rate_direct_assimilated_BC.close();

  // CFD側(時刻 t*dt_cfd)の流量
  for(size_t t = 0; t < flow_rate_cfd2.size(); ++t) {
    double time_cfd = t * dt_cfd;
    flow_rate_direct_MRI_BC << time_cfd << " " << flow_rate_cfd2[t] << "\n";
  }
  flow_rate_direct_MRI_BC.close();

  //==============================
  // 8. 相対誤差の計算 (MRI の時間点に対応)
  //==============================
  double averaged_error = 0e0;
  double averaged_error2 = 0e0;

  for(size_t i = 0; i < flow_rate_mri.size(); ++i) {
    double t_mri_sec = i * dt_mri;  // MRI 側の時刻
    // 最近傍のCFDステップ
    int idx_cfd = static_cast<int>(std::round(t_mri_sec / dt_cfd));

    // CFDのステップが範囲外ならエラー計算できないのでスキップ
    if(idx_cfd < 0 || idx_cfd >= static_cast<int>(flow_rate_cfd.size())) {
      continue;
    }
    double Qmri = flow_rate_mri[i];
    double Qcfd = flow_rate_cfd[idx_cfd];
    double Qcfd2 = flow_rate_cfd2[idx_cfd];

    // MRI側の流量が 0 のときに割り算すると発散する可能性があるのでチェック
    if(std::fabs(Qmri) < 1e-15) {
      // ゼロ割回避。必要に応じてスキップ or 別処理
      continue;
    }

    double error_percent = 100.0 * std::fabs(Qcfd - Qmri) / std::fabs(Qmri);
    double error2_percent = 100.0 * std::fabs(Qcfd2 - Qmri) / std::fabs(Qmri);

    averaged_error += error_percent;
    averaged_error2 += error2_percent;
  }

  averaged_error = averaged_error / flow_rate_mri.size();
  averaged_error2 = averaged_error2 / flow_rate_mri.size();

  std::cout << averaged_error << " " << averaged_error2 << std::endl;


  //==============================
  // 9. 正常終了
  //==============================
  return 0;
}

