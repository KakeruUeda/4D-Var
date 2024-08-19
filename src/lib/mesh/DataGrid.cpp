#include "DataGrid.h"

DataGrid::DataGrid(Config &conf, Grid &grid, SnapShot &snap)
    : grid(grid), snap(snap), vvox(conf.vvox), gauss(2), nxData(conf.nxData), nyData(conf.nyData), nzData(conf.nzData),
      lxData(conf.lxData), lyData(conf.lyData), lzData(conf.lzData), dxData(conf.dxData), dyData(conf.dyData),
      dzData(conf.dzData), nDataCellsGlobal(conf.nDataCellsGlobal), nDataNodesInCell(conf.nNodesInCellData)
{
  voxel.allocate(nzData, nyData, nxData);

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        voxel(k, j, i).v_mri.allocate(snap.nSnapShot, 3);
        voxel(k, j, i).v_cfd.allocate(snap.nSnapShot, 3);
        voxel(k, j, i).v_err.allocate(snap.nSnapShot, 3);
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

void DataGrid::initialize(Config &conf)
{
  setVoxelCenters();
  setVoxelBoundaries();
  collectCellsInVoxel();

  for(int t = 0; t < snap.nSnapShot; t++) {
    std::string file = conf.dataDir + "/data_" + std::to_string(t) + ".dat";
    importDAT(file, t);
  }
}

void DataGrid::setVoxelCenters()
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

void DataGrid::setVoxelBoundaries()
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

void DataGrid::collectCellsInVoxel()
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

void DataGrid::collectCellsInCircle(const int radious)
{
  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      if(isCellsIncludedInCircle(radious, iv, ic)) {
        voxel(iv).cells.push_back(ic);
      }
    }
  }
}

bool DataGrid::isCellsIncludedInVoxel(const int iv, const int ic)
{
  return !(voxel(iv).maxX < grid.cell(ic).minX || voxel(iv).minX > grid.cell(ic).maxX ||
           voxel(iv).maxY < grid.cell(ic).minY || voxel(iv).minY > grid.cell(ic).maxY ||
           voxel(iv).maxZ < grid.cell(ic).minZ || voxel(iv).minZ > grid.cell(ic).maxZ);
}

bool DataGrid::isCellsIncludedInCircle(double radius, int iv, int ic)
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

void DataGrid::compSmoothing()
{
  smoothing.allocate(grid.cell.nNodesInCell);
  smoothing.fillZero();

  coeff = 0.1 * std::min(std::min(dxData, dyData), dzData);

  auto sinc = [](double x) { return (x == 0) ? 1e0 : sin(PI * x) / (PI * x); };

  auto kai = [](double value, double center, double d, double a) {
    return 1 / (1 + exp(-(value - (center - 2 * d)) / a)) - 1 / (1 + exp(-(value - (center + 2 * d)) / a));
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

  for(int p = 0; p < grid.cell.nNodesInCell; p++) {
    smoothing(p) = compWeight(p, center);
  }
}

void DataGrid::weightedAverage(const int iv, const int t)
{
  if(voxel(iv).cells.size() == 0) {
    return;  // outside of the fluid domain
  }

  mt3d.nNodesInCell = grid.cell.nNodesInCell;
  mt3d.N.allocate(grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(grid.cell.nNodesInCell, 3);
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

void DataGrid::average(const int iv, const int t)
{
  if(voxel(iv).cells.size() == 0) {
    return;  // outside of the fluid domain
  }

  mt3d.nNodesInCell = grid.cell.nNodesInCell;
  mt3d.N.allocate(grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(grid.cell.nNodesInCell, 3);
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

  auto evaluateValues = [&](vector<double> &values, const int iv, const int t) {
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
          auto v_gp = mt3d.getVectorValuesGP(velCurrent);
          evaluateValues(v_gp, iv, t);
        }
      }
    }
  };

  for(int ic = 0; ic < voxel(iv).cells.size(); ic++) {
    int ivc = voxel(iv).cells[ic];
    getNodeValues(ivc, t);
    averageInVoxel(iv, t);
  }

  for(int d = 0; d < 3; d++) {
    voxel(iv).v_cfd(t, d) /= dxData * dyData * dzData;
  }
}

void DataGrid::exportDAT(const std::string &filename, const int step)
{
  std::ofstream ofs(filename);
  if(!ofs) {
    std::cerr << "Could not open file for writing: " << filename << std::endl;
    return;
  }

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        for(int d = 0; d < 3; d++) {
          ofs << voxel(k, j, i).v_cfd(step, d) << " ";
        }
        ofs << "\n";
      }
    }
  }
}

void DataGrid::importDAT(const std::string &filename, const int step)
{
  std::ifstream ifs(filename);
  if(!ifs) {
    std::cerr << "Could not open file for reading: " << filename << std::endl;
    return;
  }

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        for(int d = 0; d < 3; d++) {
          ifs >> voxel(k, j, i).v_mri(step, d);
        }
      }
    }
  }
}

void DataGrid::exportVTI(const std::string &filename, const int t)
{
  FILE *fp = fopen(filename.c_str(), "w");

  if(fp == NULL) {
    std::cout << filename << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%e %e %e\" Spacing=\"%e %e %e\">\n", 0, nxData, 0,
          nyData, 0, nzData, 0e0, 0e0, 0e0, dxData, dyData, dzData);
  fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\">\n", 0, nxData, 0, nyData, 0, nzData);
  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");
  fprintf(fp, "<CellData>\n");

  int offset = 0;
  fprintf(fp,
          "<DataArray type=\"Float32\" Name=\"vCFD\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",
          offset);
  offset += sizeof(unsigned long) + nxData * nyData * nzData * 3 * sizeof(float);
  fprintf(fp,
          "<DataArray type=\"Float32\" Name=\"vMRI\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",
          offset);
  offset += sizeof(unsigned long) + nxData * nyData * nzData * 3 * sizeof(float);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"ve\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",
          offset);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");

  fclose(fp);

  std::fstream ofs;
  ofs.open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
  unsigned long allsize = nxData * nyData * nzData * 3 * sizeof(float);

  // Write vCFD data
  float *data1 = new float[nxData * nyData * nzData * 3];
  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        int index = i + j * nxData + k * nxData * nyData;
        data1[0 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_cfd(t, 0));
        data1[1 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_cfd(t, 1));
        data1[2 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_cfd(t, 2));
      }
    }
  }

  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data1), allsize);
  delete[] data1;

  // Write vMRI data
  float *data2 = new float[nxData * nyData * nzData * 3];
  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        int index = i + j * nxData + k * nxData * nyData;
        data2[0 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_mri(t, 0));
        data2[1 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_mri(t, 1));
        data2[2 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_mri(t, 2));
      }
    }
  }

  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data2), allsize);
  delete[] data2;

  // Write ve data
  float *data3 = new float[nxData * nyData * nzData * 3];
  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        int index = i + j * nxData + k * nxData * nyData;
        data3[0 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_err(t, 0));
        data3[1 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_err(t, 1));
        data3[2 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] = static_cast<float>(voxel(k, j, i).v_err(t, 2));
      }
    }
  }

  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data3), allsize);
  delete[] data3;

  ofs.close();

  fp = fopen(filename.c_str(), "a");
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
