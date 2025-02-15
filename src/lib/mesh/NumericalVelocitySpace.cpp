#include "DataGrid.h"

void DataGrid::average_space(Array3D<double> &vt, const int iv, const int t)
{
  auto getInterpolatedVelocity = [&](const int irv, const int t, const int p) {
    int irn = refinedVoxel(irv).node[p];

    if(CFDCellId[irn] == -1) {
      return;
    }

    double ss = mt3d.xCurrent(p, 0) - grid.cell(CFDCellId[irn]).center[0];
    double tt = mt3d.xCurrent(p, 1) - grid.cell(CFDCellId[irn]).center[1];
    double uu = mt3d.xCurrent(p, 2) - grid.cell(CFDCellId[irn]).center[2];

    ss = ss / (grid.dx / 2e0);
    tt = tt / (grid.dy / 2e0);
    uu = uu / (grid.dz / 2e0);

    double eps = 1e-3;
    if(ss < -1 - eps || ss > 1 + eps) {
      throw std::runtime_error("s interpolation error.");
    } else if(tt < -1 - eps || tt > 1 + eps) {
      throw std::runtime_error("t interpolation error.");
    } else if(uu < -1 - eps || uu > 1 + eps) {
      throw std::runtime_error("u interpolation error.");
    }

    ShapeFunction3D::C3D8_N(mt3d.N, ss, tt, uu);
    for(int q = 0; q < 8; q++) {
      for(int d = 0; d < 3; d++) {
        velCurrent(p, d) += mt3d.N(q) * vt(t, grid.cell(CFDCellId[irn]).node[q], d);
      }
    }
  };

  auto getNodeValues = [&](const int irv, const int t) {
    velCurrent.fillZero();
    for(int p = 0; p < 8; p++) {
      for(int d = 0; d < 3; d++) {
        mt3d.xCurrent(p, d) = refinedVoxel(irv).x[p][d];
      }
      getInterpolatedVelocity(irv, t, p);
    }
  };

  auto evaluateValues = [&](vector<double> &values, const int iv, const int t) {
    for(int d = 0; d < 3; d++) {
      voxel(iv).v_cfd_refined(t, d) += values[d] * mt3d.vol;
    }
  };

  auto averageInVoxel = [&](const int iv, const int t) {
    for(int i1 = 0; i1 < 2; i1++) {
      for(int i2 = 0; i2 < 2; i2++) {
        for(int i3 = 0; i3 < 2; i3++) {
          mt3d.setShapesInGauss(gauss, i1, i2, i3);
          mt3d.setFactorsInGauss(gauss, i1, i2, i3);
          auto v_gp = mt3d.getVectorValuesGP(velCurrent);
          evaluateValues(v_gp, iv, t);
        }
      }
    }
  };

  for(int i = 0; i < voxel(iv).refinedVoxelId.size(); i++) {
    int irv = voxel(iv).refinedVoxelId[i];
    getNodeValues(irv, t);
    averageInVoxel(iv, t);
  }

  for(int d = 0; d < 3; d++) {
    voxel(iv).v_cfd_refined(t, d) /= dxData * dyData * dzData;
  }
}

void DataGrid::weighted_average_space(Array3D<double> &vt, const int iv, const int t)
{
  if(voxel(iv).cells.empty()) return;

  auto evaluate = [&](double &s_gp, vector<double> &v_gp, const int iv, const int t) {
    weightIntegral += s_gp * mt3d.vol;
    for(int d = 0; d < 3; d++) {
      voxel(iv).v_cfd_refined(t, d) += s_gp * v_gp[d] * mt3d.vol;
    }
  };

  auto compute_smoothing = [&](const int iv) {
    double sigma = 0.1 * std::min(std::min(dxData, dyData), dzData);

    auto sinc = [](double x) { return (x == 0) ? 1e0 : sin(PI * x) / (PI * x); };

    auto sinc_talor = [](double x) {
      if(fabs(x) < 1e-5) return 1.0 - (PI * PI * x * x) / 6.0;
      return sin(PI * x) / (PI * x);
    };

    auto kai = [](double value, double center, double d, double a) {
      return 1 / (1 + exp(-(value - (center - (d / 2e0))) / a)) - 1 / (1 + exp(-(value - (center + (d / 2e0))) / a));
    };

    auto compute_truncated_sinc = [&](int p) {
      double scale[3] = {1.0, 1.0, 1.0};
      double dx[3] = {dxData, dyData, dzData};
      double lx[3] = {lxData, lyData, lzData};

      for(int i = 0; i < 3; i++) {
        if(mt3d.xCurrent(p, i) < dx[i] / 2.0 || mt3d.xCurrent(p, i) > lx[i] - dx[i] / 2.0) {
          scale[i] = 4.0;
        }
      }

      double sinc_x = sinc((scale[0] * (mt3d.xCurrent(p, 0) - voxel(iv).center[0])) / dxData);
      double sinc_y = sinc((scale[1] * (mt3d.xCurrent(p, 1) - voxel(iv).center[1])) / dyData);
      double sinc_z = sinc((scale[2] * (mt3d.xCurrent(p, 2) - voxel(iv).center[2])) / dzData);

      double kai_x = kai(scale[0] * mt3d.xCurrent(p, 0), scale[0] * voxel(iv).center[0], 4 * dxData, sigma);
      double kai_y = kai(scale[1] * mt3d.xCurrent(p, 1), scale[1] * voxel(iv).center[1], 4 * dyData, sigma);
      double kai_z = kai(scale[2] * mt3d.xCurrent(p, 2), scale[2] * voxel(iv).center[2], 4 * dzData, sigma);

      return sinc_x * sinc_y * sinc_z * kai_x * kai_y * kai_z;
    };

    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      smoothing(p) = compute_truncated_sinc(p);
    }
  };

  weightIntegral = 0e0;

  for(int ic = 0; ic < voxel(iv).cells.size(); ic++) {
    int ivc = voxel(iv).cells[ic];

    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      for(int d = 0; d < 3; d++) {
        velCurrent(p, d) = vt(t, grid.cell(ivc).node[p], d);
        mt3d.xCurrent(p, d) = grid.cell(ivc).x[p][d];
      }
    }

    compute_smoothing(iv);

    for(int i1 = 0; i1 < 2; i1++) {
      for(int i2 = 0; i2 < 2; i2++) {
        for(int i3 = 0; i3 < 2; i3++) {
          mt3d.setShapesInGauss(gauss, i1, i2, i3);
          mt3d.setFactorsInGauss(gauss, i1, i2, i3);

          auto s_gp = mt3d.getScalarValueGP(smoothing);
          auto v_gp = mt3d.getVectorValuesGP(velCurrent);

          weightIntegral += s_gp * mt3d.vol;
          for(int d = 0; d < 3; d++) {
            voxel(iv).v_cfd_refined(t, d) += s_gp * v_gp[d] * mt3d.vol;
          }
        }
      }
    }
  }

  if(weightIntegral == 0) {
    throw std::runtime_error("Weight integral is zero");
  }

  for(int d = 0; d < 3; d++) {
    voxel(iv).v_cfd_refined(t, d) /= weightIntegral;
  }
}

void DataGrid::linear_interpolate_space(Array3D<double> &vt, const int iv, const int t)
{
  if(voxel(iv).CFDCellId == -1) return;

  mt3d.nNodesInCell = grid.cell.nNodesInCell;
  mt3d.N.allocate(grid.cell.nNodesInCell);
  velCurrent.allocate(grid.cell.nNodesInCell, 3);

  int ic = voxel(iv).CFDCellId;

  double ss = voxel(iv).center[0] - grid.cell(ic).center[0];
  double tt = voxel(iv).center[1] - grid.cell(ic).center[1];
  double uu = voxel(iv).center[2] - grid.cell(ic).center[2];

  ss = ss / (grid.dx / 2e0);
  tt = tt / (grid.dy / 2e0);
  uu = uu / (grid.dz / 2e0);

  double eps = 1e-3;
  if(ss < -1 - eps || ss > 1 + eps) {
    throw std::runtime_error("s interpolation error.");
  } else if(tt < -1 - eps || tt > 1 + eps) {
    throw std::runtime_error("t interpolation error.");
  } else if(uu < -1 - eps || uu > 1 + eps) {
    throw std::runtime_error("u interpolation error.");
  }

  ShapeFunction3D::C3D8_N(mt3d.N, ss, tt, uu);

  for(int p = 0; p < mt3d.nNodesInCell; p++) {
    for(int d = 0; d < 3; d++) {
      velCurrent(p, d) = vt(t, grid.cell(ic).node[p], d);
    }
  }

  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < mt3d.nNodesInCell; p++) {
      voxel(iv).v_cfd_refined(t, d) += mt3d.N(p) * velCurrent(p, d);
    }
  }
}