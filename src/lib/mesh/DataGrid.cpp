#include "DataGrid.h"

DataGrid::DataGrid(Config &conf, Grid &grid, SnapShot &snap)
    : grid(grid), snap(snap), vel_space(conf.vel_space), vel_time(conf.vel_time), gauss(2), g1(1), nxData(conf.nxData),
      nyData(conf.nyData), nzData(conf.nzData), lxData(conf.lxData), lyData(conf.lyData), lzData(conf.lzData),
      dxData(conf.dxData), dyData(conf.dyData), dzData(conf.dzData), nDataCellsGlobal(conf.nDataCellsGlobal),
      nDataNodesInCell(conf.nNodesInCellData), dt_mri(conf.dt_mri), dt_cfd(conf.dt), n_mri_step(conf.nSnapShot),
      n_cfd_step(conf.timeMax)
{
  voxel.allocate(nzData, nyData, nxData);

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        voxel(k, j, i).v_mri.allocate(snap.nSnapShot, 3);
        voxel(k, j, i).v_cfd.allocate(snap.nSnapShot, 3);
        voxel(k, j, i).v_err.allocate(snap.nSnapShot, 3);
        voxel(k, j, i).v_mri_refined.allocate(snap.timeMax, 3);
        voxel(k, j, i).v_cfd_refined.allocate(snap.timeMax, 3);
        voxel(k, j, i).v_err_refined.allocate(snap.timeMax, 3);
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
  for(int t = 0; t < snap.nSnapShot; t++) {
    std::string file = conf.dataDir + "/data_" + std::to_string(t) + ".dat";
    importDAT(file, t);
  }

  std::string file_mask = conf.dataDir + "/mask.dat";
  importMask(file_mask);

  for(int iv = 0; iv < nDataCellsGlobal; ++iv) {
    voxel(iv).subId = conf.voxelId[iv];
  }

  local_voxel_info.resize(n_cfd_step * nDataCellsGlobal * 3, 0e0);
  global_voxel_info.resize(local_voxel_info.size(), 0e0);

  grid.cell.getBoundaries();
  grid.cell.getCenterCoordinates();
  setVoxelCenters();
  setBoundingBox();
  //setVoxelBoundaries();
  //collectCellsInVoxel();

  setRefinedGrid();
  collectRefinedVoxelId();
  collectCFDCellId();
  //collectVoxelCenterCFDCell();

  collectCellsInCircle(2 * dxData);
  collectCFDTimeStep();

  extract_slice();

  mt1d.nNodesInCell = 2;
  mt1d.N.allocate(2);
  mt1d.dNdr.allocate(2);
  mt1d.xCurrent.allocate(2);
  vel_current_time.allocate(2, 3);
  smoothing_time.allocate(2);
  smoothing_time.fillZero();

  mt3d.nNodesInCell = grid.cell.nNodesInCell;
  mt3d.N.allocate(grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(grid.cell.nNodesInCell, 3);
  mt3d.xCurrent.allocate(grid.cell.nNodesInCell, 3);
  velCurrent.allocate(grid.cell.nNodesInCell, 3);
  smoothing.allocate(grid.cell.nNodesInCell);
  smoothing.fillZero();
}

void DataGrid::collectCFDTimeStep()
{
  double time_min, time_max;

  cfd_time_steps_id.resize(n_mri_step);

  for(int t1 = 0; t1 < n_mri_step; t1++) {
    double time_mri = t1 * dt_mri;
    time_min = time_mri - dt_mri / 2e0;
    time_max = time_mri + dt_mri / 2e0;

    if(t1 == 0) time_min = -1e10;
    if(t1 == n_mri_step - 1) time_max = 1e10;

    for(int t2 = 0; t2 < n_cfd_step; t2++) {
      double time_cfd = t2 * dt_cfd;

      if(time_cfd >= time_min && time_cfd <= time_max) {
        cfd_time_steps_id[t1].push_back(t2);
      }
    }
  }

  // if(mpi.myId == 0) {
  //   ofstream ofs("id.dat");
  //   for(int t1 = 0; t1 < n_mri_step; t1++) {
  //     for(int t2 = 0; t2 < cfd_time_steps_id[t1].size(); t2++) {
  //       ofs << cfd_time_steps_id[t1][t2] << " ";
  //     }
  //     ofs << std::endl;
  //   }
  //   ofs.close();
  // }
}

void DataGrid::setRefinedGrid()
{
  refinementFactor = 2;

  nxRefinedData = nxData * refinementFactor;
  nyRefinedData = nyData * refinementFactor;
  nzRefinedData = nzData * refinementFactor;

  dxRefinedData = dxData / refinementFactor;
  dyRefinedData = dyData / refinementFactor;
  dzRefinedData = dzData / refinementFactor;

  nRefinedCellGlobal = nxRefinedData * nyRefinedData * nzRefinedData;
  nRefinedNodeGlobal = (nxRefinedData + 1) * (nyRefinedData + 1) * (nzRefinedData + 1);

  refinedVoxel.allocate(nzRefinedData, nyRefinedData, nxRefinedData);

  for(int k = 0; k < nzRefinedData; k++) {
    for(int j = 0; j < nyRefinedData; j++) {
      for(int i = 0; i < nxRefinedData; i++) {
        refinedVoxel(k, j, i).node.resize(8, 0);
        for(int p = 0; p < 8; p++) {
          refinedVoxel(k, j, i).node[p] = Config::setStrNode(i, j, k, p, nxRefinedData, nyRefinedData, nzRefinedData);
        }
      }
    }
  }

  refinedCoord.resize(nRefinedNodeGlobal, std::vector<double>(3, 0e0));

  for(int k = 0; k < nzRefinedData + 1; k++) {
    for(int j = 0; j < nyRefinedData + 1; j++) {
      for(int i = 0; i < nxRefinedData + 1; i++) {
        int in = k * (nyRefinedData + 1) * (nxRefinedData + 1) + j * (nxRefinedData + 1) + i;
        for(int d = 0; d < 3; d++) {
          refinedCoord[in][d] = Config::setStrCoordinate(i, j, k, d, dxRefinedData, dyRefinedData, dzRefinedData);
        }
      }
    }
  }

  for(int k = 0; k < nzRefinedData + 1; k++) {
    for(int j = 0; j < nyRefinedData + 1; j++) {
      for(int i = 0; i < nxRefinedData + 1; i++) {
        int in = k * (nyRefinedData + 1) * (nxRefinedData + 1) + j * (nxRefinedData + 1) + i;
        refinedCoord[in][0] += xOrigin;
        refinedCoord[in][1] += yOrigin;
        refinedCoord[in][2] += zOrigin;
      }
    }
  }

  for(int irv = 0; irv < nRefinedCellGlobal; irv++) {
    refinedVoxel(irv).x.resize(8, std::vector<double>(3, 0e0));
    for(int p = 0; p < 8; p++) {
      int in = refinedVoxel(irv).node[p];
      for(int d = 0; d < 3; d++) {
        refinedVoxel(irv).x[p][d] = refinedCoord[in][d];
      }
    }
  }
}

void DataGrid::collectRefinedVoxelId()
{
  auto collectId = [&](int k, int j, int i) {
    for(int u = k * refinementFactor; u < (k + 1) * refinementFactor; u++) {
      for(int t = j * refinementFactor; t < (j + 1) * refinementFactor; t++) {
        for(int s = i * refinementFactor; s < (i + 1) * refinementFactor; s++) {
          int irv = s + t * nxRefinedData + u * nxRefinedData * nyRefinedData;
          voxel(k, j, i).refinedVoxelId.push_back(irv);
        }
      }
    }
  };

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        collectId(k, j, i);
      }
    }
  }
}

void DataGrid::setBoundingBox()
{
  box.minX = grid.cell(0).minX;
  box.maxX = grid.cell(0).maxX;
  box.minY = grid.cell(0).minY;
  box.maxY = grid.cell(0).maxY;
  box.minZ = grid.cell(0).minZ;
  box.maxZ = grid.cell(0).maxZ;

  for(int ic = 1; ic < grid.cell.nCellsGlobal; ic++) {
    if(grid.cell(ic).minX < box.minX) {
      box.minX = grid.cell(ic).minX;
    }
    if(grid.cell(ic).maxX > box.maxX) {
      box.maxX = grid.cell(ic).maxX;
    }

    if(grid.cell(ic).minY < box.minY) {
      box.minY = grid.cell(ic).minY;
    }
    if(grid.cell(ic).maxY > box.maxY) {
      box.maxY = grid.cell(ic).maxY;
    }

    if(grid.cell(ic).minZ < box.minZ) {
      box.minZ = grid.cell(ic).minZ;
    }
    if(grid.cell(ic).maxZ > box.maxZ) {
      box.maxZ = grid.cell(ic).maxZ;
    }
  }
}

bool DataGrid::isRefinedNodeOutside(const int in)
{
  return (refinedCoord[in][0] < box.minX || refinedCoord[in][0] > box.maxX || refinedCoord[in][1] < box.minY ||
          refinedCoord[in][1] > box.maxY || refinedCoord[in][2] < box.minZ || refinedCoord[in][2] > box.maxZ);
}

void DataGrid::collectCFDCellId()
{
  CFDCellId.resize(nRefinedNodeGlobal, -1);

  Octree octree(grid.cell, xOrigin, xOrigin + lxData, yOrigin, yOrigin + lyData, zOrigin, zOrigin + lzData, 3);

  for(int in = 0; in < nRefinedNodeGlobal; in++) {
    bool flag = false;

    std::vector<int> candidateCells = octree.getCandidateCells(refinedCoord[in]);

    for(int ic : candidateCells) {
      if(isRefinedNodeIncludedInCFDCell(in, ic)) {
        CFDCellId[in] = ic;
        flag = true;
      }
      if(flag) {
        break;
      }
    }
  }
  // if(mpi.myId == 0) {
  //   std::string file = "id.vti";
  //   EXPORT::exportScalarPointDataVTI(file, "id", CFDCellId, nxRefinedData, nyRefinedData, nzRefinedData, dxRefinedData, dyRefinedData, dzRefinedData);
  // }
}

void DataGrid::collectVoxelCenterCFDCell()
{
  Octree octree(grid.cell, 0e0, lxData, 0e0, lyData, 0e0, lzData, 3);

  std::vector<std::vector<double>> coord(nDataCellsGlobal, std::vector<double>(3, 0e0));

  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    bool flag = false;

    for(int d = 0; d < 3; d++) {
      coord[iv][d] = voxel(iv).center[d];
    }

    std::vector<int> candidateCells = octree.getCandidateCells(coord[iv]);

    voxel(iv).CFDCellId = -1;

    for(int ic : candidateCells) {
      if(isVoxelCenterIncludedInCFDCell(iv, ic)) {
        voxel(iv).CFDCellId = ic;
        flag = true;
      }
      if(flag) {
        break;
      }
    }
  }

  // std::vector<int> id;

  // for(int iv = 0; iv < nDataCellsGlobal; iv++) {
  //   id.push_back(voxel(iv).CFDCellId);
  // }
  // if(mpi.myId == 0) {
  //   std::string file = "id.vti";
  //   EXPORT::exportScalarCellDataVTI(file, "id", id, nxData, nyData, nzData, dxData, dyData, dzData);
  // }
}

bool DataGrid::isRefinedNodeIncludedInCFDCell(const int in, const int ic)
{
  return !(refinedCoord[in][0] < grid.cell(ic).minX || refinedCoord[in][0] > grid.cell(ic).maxX ||
           refinedCoord[in][1] < grid.cell(ic).minY || refinedCoord[in][1] > grid.cell(ic).maxY ||
           refinedCoord[in][2] < grid.cell(ic).minZ || refinedCoord[in][2] > grid.cell(ic).maxZ);
}

bool DataGrid::isVoxelCenterIncludedInCFDCell(const int iv, const int ic)
{
  return !(voxel(iv).center[0] < grid.cell(ic).minX || voxel(iv).center[0] > grid.cell(ic).maxX ||
           voxel(iv).center[1] < grid.cell(ic).minY || voxel(iv).center[1] > grid.cell(ic).maxY ||
           voxel(iv).center[2] < grid.cell(ic).minZ || voxel(iv).center[2] > grid.cell(ic).maxZ);
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
  //grid.cell.getBoundaries();
  grid.cell.getCenterCoordinates();

  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      if(isCellCenterIncludedInVoxel(iv, ic)) {
        voxel(iv).cells.push_back(ic);
      }
    }
  }
}

void DataGrid::collectCellsInCircle(double radious)
{
  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    if(voxel(iv).mask < 0.5) continue;
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      if(isCellIncludedInCircle(radious, iv, ic)) {
        voxel(iv).cells.push_back(ic);
      }
    }
  }
}

bool DataGrid::isCellCenterIncludedInVoxel(const int iv, const int ic)
{
  return !(voxel(iv).maxX < grid.cell(ic).center[0] || voxel(iv).minX > grid.cell(ic).center[0] ||
           voxel(iv).maxY < grid.cell(ic).center[1] || voxel(iv).minY > grid.cell(ic).center[1] ||
           voxel(iv).maxZ < grid.cell(ic).center[2] || voxel(iv).minZ > grid.cell(ic).center[2]);
}

bool DataGrid::isCellBoundaryIncludedInVoxel(const int iv, const int ic)
{
  return !(voxel(iv).maxX < grid.cell(ic).minX || voxel(iv).minX > grid.cell(ic).maxX ||
           voxel(iv).maxY < grid.cell(ic).minY || voxel(iv).minY > grid.cell(ic).maxY ||
           voxel(iv).maxZ < grid.cell(ic).minZ || voxel(iv).minZ > grid.cell(ic).maxZ);
}

bool DataGrid::isCellIncludedInCircle(double radius, int iv, int ic)
{
  double radius2 = radius * radius;

  for(int p = 0; p < grid.cell.nNodesInCell; p++) {
    double distance = 0e0;
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

/**
 * @brief Compute continuous v_mri
 */
void DataGrid::extract_slice()
{
  for(int k = 0; k < nzData; ++k) {
    for(int i = 0; i < nxData; ++i) {
      int iv = (k * nxData * nyData) + i;
      if(voxel(iv).mask < 1e-12) continue;
      x_slice.push_back(i * dxData + dxData / 2);
      z_slice.push_back(k * dzData + dzData / 2);
    }
  }

  u_mri_slice.resize(n_mri_step);
  v_mri_slice.resize(n_mri_step);
  w_mri_slice.resize(n_mri_step);

  for(int t = 0; t < n_mri_step; t++) {
    for(int k = 0; k < nzData; ++k) {
      for(int i = 0; i < nxData; ++i) {
        int iv = (k * nxData * nyData) + i;
        if(voxel(iv).mask < 1e-12) continue;
        if(voxel(iv).mask < 0.9999999) {
          u_mri_slice[t].push_back(0e0);
          v_mri_slice[t].push_back(0e0);
          w_mri_slice[t].push_back(0e0);
        } else {
          u_mri_slice[t].push_back(voxel(iv).v_mri(t, 0));
          v_mri_slice[t].push_back(voxel(iv).v_mri(t, 1));
          w_mri_slice[t].push_back(voxel(iv).v_mri(t, 2));
        }
      }
    }
  }
}

/**
 * @brief Compute continuous v_mri
 */
void DataGrid::comp_v_mri_refined(const double dt, const int timeMax)
{
  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    std::vector<double> x, y1, y2, y3;

    for(int t = 0; t < snap.nSnapShot; t++) {
      double dp = t * dt * snap.snapInterval;
      x.push_back(dp);
      y1.push_back(voxel(iv).v_mri(t, 0));
      y2.push_back(voxel(iv).v_mri(t, 1));
      y3.push_back(voxel(iv).v_mri(t, 2));
    }

    vector<Spline::Coefficients> cf_x = Spline::compCoefficients(x, y1);
    vector<Spline::Coefficients> cf_y = Spline::compCoefficients(x, y2);
    vector<Spline::Coefficients> cf_z = Spline::compCoefficients(x, y3);

    for(int t = 0; t < timeMax; t++) {
      double p = t * dt;
      voxel(iv).v_mri_refined(t, 0) = Spline::evaluate(cf_x, p);
      voxel(iv).v_mri_refined(t, 1) = Spline::evaluate(cf_y, p);
      voxel(iv).v_mri_refined(t, 2) = Spline::evaluate(cf_z, p);
    }
  }
}

/**
 * @brief Compute continuous v_err
 */
void DataGrid::comp_v_err_refined()
{
  for(int iv = 0; iv < nDataCellsGlobal; iv++) {
    std::vector<double> x, y1, y2, y3;

    for(int t = 0; t < n_mri_step; t++) {
      double p = t * dt_mri;
      x.push_back(p);
      y1.push_back(voxel(iv).v_err(t, 0));
      y2.push_back(voxel(iv).v_err(t, 1));
      y3.push_back(voxel(iv).v_err(t, 2));
    }

    vector<Spline::Coefficients> cf_x = Spline::compCoefficients(x, y1);
    vector<Spline::Coefficients> cf_y = Spline::compCoefficients(x, y2);
    vector<Spline::Coefficients> cf_z = Spline::compCoefficients(x, y3);

    for(int t = 0; t < n_cfd_step; t++) {
      double p = t * dt_cfd;
      voxel(iv).v_err_refined(t, 0) = Spline::evaluate(cf_x, p);
      voxel(iv).v_err_refined(t, 1) = Spline::evaluate(cf_y, p);
      voxel(iv).v_err_refined(t, 2) = Spline::evaluate(cf_z, p);
    }
  }
}

void DataGrid::gather_voxel_info()
{
  int index = 0;
  for(int t = 0; t < n_cfd_step; t++) {
    for(int iv = 0; iv < nDataCellsGlobal; iv++) {
      for(int d = 0; d < 3; d++) {
        local_voxel_info[index] = voxel(iv).v_cfd_refined(t, d);
        index++;
      }
    }
  }

  std::fill(global_voxel_info.begin(), global_voxel_info.end(), 0e0);

  MPI_Allreduce(local_voxel_info.data(),   // send buffer
                global_voxel_info.data(),  // receive buffer
                local_voxel_info.size(),   // n_data
                MPI_DOUBLE,                // data type
                MPI_SUM,                   // summation
                MPI_COMM_WORLD             // communicator
  );

  index = 0;
  for(int t = 0; t < n_cfd_step; t++) {
    for(int iv = 0; iv < nDataCellsGlobal; iv++) {
      for(int d = 0; d < 3; d++) {
        voxel(iv).v_cfd_refined(t, d) = global_voxel_info[index];
        index++;
      }
    }
  }
}

void DataGrid::exportVelCFD_DAT(const std::string &filename, const int step)
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

void DataGrid::exportVelMRI_DAT(const std::string &filename, const int step)
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
          ofs << voxel(k, j, i).v_mri(step, d) << " ";
        }
        ofs << "\n";
      }
    }
  }
}

void DataGrid::exportVelError_DAT(const std::string &filename, const int step)
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
          ofs << voxel(k, j, i).v_err(step, d) << " ";
        }
        ofs << "\n";
      }
    }
  }
}

void DataGrid::exportVelCFD_t_DAT(const std::string &filename, const int step)
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
          ofs << voxel(k, j, i).v_cfd_refined(step, d) << " ";
        }
        ofs << "\n";
      }
    }
  }
}

void DataGrid::exportVelMRI_t_DAT(const std::string &filename, const int step)
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
          ofs << voxel(k, j, i).v_mri_refined(step, d) << " ";
        }
        ofs << "\n";
      }
    }
  }
}

void DataGrid::exportVelError_t_DAT(const std::string &filename, const int step)
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
          ofs << voxel(k, j, i).v_err_refined(step, d) << " ";
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

void DataGrid::importMask(const std::string &filename)
{
  std::ifstream ifs(filename);
  if(!ifs) {
    std::cerr << "Could not open file for reading: " << filename << std::endl;
    return;
  }

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        ifs >> voxel(k, j, i).mask;
      }
    }
  }
}

void DataGrid::exportMaskVTI(const std::string &filename)
{
  FILE *fp;
  fp = fopen(filename.c_str(), "w");

  if(fp == NULL) {
    std::cerr << filename << " open error" << std::endl;
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
  fprintf(fp,
          "<DataArray type=\"Float32\" Name=\"%mask\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n");
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  unsigned long allsize = nxData * nyData * nzData * sizeof(float);

  float *data = new float[nxData * nyData * nzData];

  for(int k = 0; k < nzData; k++) {
    for(int j = 0; j < nyData; j++) {
      for(int i = 0; i < nxData; i++) {
        int iv = i + j * nxData + k * nxData * nyData;
        data[iv] = static_cast<float>(voxel(k, j, i).mask);
      }
    }
  }

  std::fstream ofs;
  ofs.open(filename.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data), allsize);
  ofs.close();

  delete[] data;

  fp = fopen(filename.c_str(), "a");
  fprintf(fp, "\n");
  fprintf(fp, "</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
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

void DataGrid::exportVTI_t(const std::string &filename, const int t)
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
        data1[0 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_cfd_refined(t, 0));
        data1[1 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_cfd_refined(t, 1));
        data1[2 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_cfd_refined(t, 2));
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
        data2[0 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_mri_refined(t, 0));
        data2[1 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_mri_refined(t, 1));
        data2[2 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_mri_refined(t, 2));
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
        data3[0 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_err_refined(t, 0));
        data3[1 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_err_refined(t, 1));
        data3[2 + i * 3 + j * nxData * 3 + k * nxData * nyData * 3] =
            static_cast<float>(voxel(k, j, i).v_err_refined(t, 2));
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