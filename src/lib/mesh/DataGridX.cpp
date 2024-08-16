#include "DataGridX.h"

DataGridX::DataGridX(Config &conf, Grid &grid, SnapShot &snap)
    : grid(grid), snap(snap), vvox(conf.vvox), gauss(2), nxData(conf.nxData), nyData(conf.nyData), nzData(conf.nzData),
      lxData(conf.lxData), lyData(conf.lyData), lzData(conf.lzData), dxData(conf.dxData), dyData(conf.dyData),
      dzData(conf.dzData), nDataCellsGlobal(conf.nDataCellsGlobal), nDataNodesInCell(conf.nNodesInCellData)
{
  voxel.allocate(nzData, nyData, nxData);
  for(int z = 0; z < nzData; ++z) {
    for(int y = 0; y < nyData; ++y) {
      for(int x = 0; x < nxData; ++x) {
        voxel(z, y, x).v_mri.allocate(snap.nSnapShot, 3);
        voxel(z, y, x).v_cfd.allocate(snap.nSnapShot, 3);
        voxel(z, y, x).v_err.allocate(snap.nSnapShot, 3);
      }
    }
  }
  if(conf.app == Application::VOXELDATACREATION) {
    xOrigin = conf.xOrigin;
    yOrigin = conf.yOrigin;
    zOrigin = conf.zOrigin;
  } else if(conf.app == Application::FDVAR) {
    xOrigin = 0e0;
    yOrigin = 0e0;
    zOrigin = 0e0;
  }
}

void DataGridX::setVoxelCenters()
{
  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        voxel(k, j, i).center[0] = xOrigin + (5e-1 + i) * dxData;
        voxel(k, j, i).center[1] = yOrigin + (5e-1 + j) * dyData;
        voxel(k, j, i).center[2] = zOrigin + (5e-1 + k) * dzData;
      }
    }
  }
}

void DataGridX::setVoxelBoundaries()
{
  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        voxel(k, j, i).minX = voxel(k, j, i).center[0] - 5e-1 * dxData;
        voxel(k, j, i).minY = voxel(k, j, i).center[1] - 5e-1 * dyData;
        voxel(k, j, i).minZ = voxel(k, j, i).center[2] - 5e-1 * dzData;
        voxel(k, j, i).maxX = voxel(k, j, i).center[0] + 5e-1 * dxData;
        voxel(k, j, i).maxY = voxel(k, j, i).center[1] + 5e-1 * dyData;
        voxel(k, j, i).maxZ = voxel(k, j, i).center[2] + 5e-1 * dzData;
      }
    }
  }
}

void DataGridX::collectCellsInVoxel()
{
  grid.cell.getBoundaries();
  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      if(isCellsIncludedInVoxel(iv, ic)) {
        voxel(iv).cells.push_back(ic);
      }
    }
  }
}

void DataGridX::collectCellsInCircle(const int radious)
{
  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      if(isCellsIncludedInCircle(radious, iv, ic)) {
        voxel(iv).cells.push_back(ic);
      }
    }
  }
}

bool DataGridX::isCellsIncludedInVoxel(const int iv, const int ic)
{
  return !(voxel(iv).maxX < grid.cell(ic).minX || voxel(iv).minX > grid.cell(ic).maxX ||
           voxel(iv).maxY < grid.cell(ic).minY || voxel(iv).minY > grid.cell(ic).maxY ||
           voxel(iv).maxZ < grid.cell(ic).minZ || voxel(iv).minZ > grid.cell(ic).maxZ);
}

bool DataGridX::isCellsIncludedInCircle(double radius, int iv, int ic)
{
  double radius2 = radius * radius;

  for(int p = 0; p < grid.cell.nNodesInCell; p++) {
    double distance = 0.0;
    for(int d = 0; d < 3; ++d) {
      double diff = grid.cell(ic).x[p][d] - voxel(iv).center[d];
      distance += diff * diff;
    }
    if(distance < radius2) {
      return true;
    }
  }
  return false;
}

void DataGridX::compSmoothing()
{
  smoothing.allocate(grid.cell.nNodesInCell);
  smoothing.fillZero();

  coeff = 0.1 * std::min(std::min(dxData, dyData), dzData);

  auto sinc = [](double x) { return (x == 0) ? 1e0 : sin(PI * x) / (PI * x); };

  auto kai = [](double value, double center, double d, double a) {
    return 1 / (1 + exp(-(value - (center - 2 * d)) / a)) - 1 
             / (1 + exp(-(value - (center + 2 * d)) / a));
  };

  auto compCenter = [&]() {
    double center[3] = {0e0, 0e0, 0e0};
    ShapeFunction3D::C3D8_N(mt3d.N, 0e0, 0e0, 0e0);
    for(int d = 0; d < 3; d++) {
      for(int p = 0; p < grid.cell.nNodesInCell; p++) {
        center[d] += mt3d.N(p) * mt3d.xCurrent(p, d);
      }
    }
    return std::array<double, 3>{center[0], center[1], center[2]};
  };

  auto compWeight = [&](int p, const std::array<double, 3> &center) {
    double sinc_x = sinc((mt3d.xCurrent(p, 0) - center[0]) / grid.dx);
    double sinc_y = sinc((mt3d.xCurrent(p, 0) - center[1]) / grid.dy);
    double sinc_z = sinc((mt3d.xCurrent(p, 0) - center[2]) / grid.dz);

    double kai_x = kai(mt3d.xCurrent(p, 0), center[0], grid.dx, coeff);
    double kai_y = kai(mt3d.xCurrent(p, 0), center[1], grid.dy, coeff);
    double kai_z = kai(mt3d.xCurrent(p, 0), center[2], grid.dz, coeff);

    return sinc_x * sinc_y * sinc_z * kai_x * kai_y * kai_z;
  };

  auto center = compCenter();

  if(vvox == VoxelVelocity::POINTSPREAD) {
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      smoothing(p) = compWeight(p, center);
    }
  } else if(vvox == VoxelVelocity::AVERAGE) {
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      smoothing(p) = 1.0;
    }
  }
}

void DataGridX::average(const int iv, const int t)
{
  mt3d.nNodesInCell = grid.cell.nNodesInCell;
  mt3d.N.allocate(grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(grid.cell.nNodesInCell, 3);
  mt3d.xCurrent.allocate(grid.cell.nNodesInCell, 3);
  velCurrent.allocate(grid.cell.nNodesInCell, 3);

  auto getNodeValues = [&](const int ivc, const int t) {
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      for(int d = 0; d < 3; d++) {
        velCurrent(p, d) = snap.vSnap(t, grid.cell(ivc).node[p], d);
        mt3d.xCurrent(p, d) = grid.cell(ivc).x[p][d];
      }
    }
  };

  auto evaluateValues = [&](double &value, vector<double> &values, const int iv, const int t) {
    weightIntegral += value * mt3d.vol;
    for(int d = 0; d < 3; d++) {
      voxel(iv).v_cfd(t, d) += values[d] * mt3d.vol;
    }
  };

  auto averageInVoxel = [&](const int iv, const int t) {
    for(int i1 = 0; i1 < 2; i1++) {
      for(int i2 = 0; i2 < 2; i2++) {
        for(int i3 = 0; i3 < 2; i3++) {
          mt3d.setShapesInGauss(gauss, i1, i2, i3);
          mt3d.setFactorsInGauss(gauss, i1, i2, i3);
          auto s_gp = mt3d.getScalarValueGP(smoothing);
          auto v_gp = mt3d.getVectorValuesGP(velCurrent);
          evaluateValues(s_gp, v_gp, iv, t);
        }
      }
    }
  };

  weightIntegral = 0e0;

  for(int ic = 0; ic < voxel(iv).cells.size(); ic++) {
    int ivc = voxel(iv).cells[ic];
    compSmoothing();
    getNodeValues(ivc, t);
    averageInVoxel(iv, t);
  }

  if(weightIntegral == 0) {
    throw std::runtime_error("Weight integral is zero");
  }

  for(int d = 0; d < 3; d++) {
    voxel(iv).v_cfd(t, d) /= weightIntegral;
  }
}