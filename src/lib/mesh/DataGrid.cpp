/**
 * @file DataGrid.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DataGrid.h"

DataGrid::DataGrid(Config &conf) :
nx(conf.nxData), ny(conf.nyData), nz(conf.nzData),
lx(conf.lxData), ly(conf.lyData), lz(conf.lzData),
dx(conf.dxData), dy(conf.dyData), dz(conf.dzData),
nData(conf.nData), nSnapShot(conf.nSnapShot), dim(conf.dim),
snapInterval(conf.snapInterval),
nCellsGlobal(conf.nCellsDataGlobal),
nNodesInCell(conf.nNodesInCellData),
data(conf.nCellsDataGlobal)
{
  for (int ic = 0; ic < conf.nCellsDataGlobal; ic++)
  {
    data[ic].dx = dx;
    data[ic].dy = dy;
    data[ic].dz = dz;
    data[ic].nNodesInCell = nNodesInCell;
  }
  if (conf.app == Application::VOXELDATA)
  {
    xOrigin = conf.xOrigin;
    yOrigin = conf.yOrigin;
    zOrigin = conf.zOrigin;
    for (int ic = 0; ic < conf.nCellsDataGlobal; ic++)
    {
      VecTool::resize(data[ic].vCFD, conf.nSnapShot, conf.dim);
      VecTool::resize(data[ic].vMRI, conf.nSnapShot, conf.dim);
      VecTool::resize(data[ic].ve, conf.nSnapShot, conf.dim);
      VecTool::resize(data[ic].center, conf.dim);
    }
  }
  else if (conf.app == Application::FDVAR)
  {
    for (int ic = 0; ic < conf.nCellsDataGlobal; ic++)
    {
      VecTool::resize(data[ic].vCFD, conf.nSnapShot, conf.dim);
      VecTool::resize(data[ic].vMRI, conf.nSnapShot, conf.dim);
      VecTool::resize(data[ic].ve, conf.nSnapShot, conf.dim);
      VecTool::resize(data[ic].center, conf.dim);
    }
    VecTool::resize(vEX, nz + 2, ny + 2, nx + 2, conf.dim);
  }
}

void DataGrid::initialize(Config &conf, Node &node, Cell &cell, const int &dim)
{
  range = 5e-1 * sqrt(dx * dx + dy * dy + dz * dz);

  for (int t = 0; t < conf.nSnapShot; t++)
  {
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int i = 0; i < nx; i++)
        {
          for (int d = 0; d < dim; d++)
          {
            data[k * nx * ny + j * nx + i].vMRI[t][d] 
            = conf.velocityData[t][k * nx * ny + j * nx + i][d];
          }
        }
      }
    }
  }

  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        data[k * nx * ny + j * nx + i].center[0] = (5e-1 + i) * dx;
        data[k * nx * ny + j * nx + i].center[1] = (5e-1 + j) * dy;
        data[k * nx * ny + j * nx + i].center[2] = (5e-1 + k) * dz;
        data[k * nx * ny + j * nx + i].setNearCell(node, cell, range, dim);
      }
    }
  }

  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        data[k * nx * ny + j * nx + i].setCellOnCenterPoint(node, cell, dim);
      }
    }
  }
}

void DataGrid::compEdgeValue(const int t)
{
  for (int k = 0; k < nz + 2; k++)
  {
    for (int j = 0; j < ny + 2; j++)
    {
      for (int i = 0; i < nx + 2; i++)
      {
        for (int d = 0; d < dim; d++)
        {
          vEX[k][j][i][d] = 0e0;
        }
      }
    }
  }
  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        for (int d = 0; d < dim; d++)
        {
          vEX[k + 1][j + 1][i + 1][d] 
          = data[k * nx * ny + j * nx + i].ve[t][d];
        }
      }
    }
  }
}

void VoxelInfo::setNearCell(Node &node, Cell &cell, const double range, const int dim)
{
  double distance;
  double diff[dim];
  bool flag;

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    flag = false;
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      distance = 0e0;
      for (int d = 0; d < dim; d++)
      {
        diff[d] = node.x[cell(ic).node[p]][d] - center[d];
        distance += diff[d] * diff[d];
      }
      distance = sqrt(distance);
      if (distance < range)
        flag = true;
    }
    if (flag)
      cellChildren.push_back(ic);
  }
}

void VoxelInfo::setCellOnCenterPoint(Node &node, Cell &cell, const int dim)
{
  double distance;
  std::vector<double> diff(dim);
  std::vector<double> minDiff(dim);
  std::vector<double> point(dim);

  double dx, dy, dz;

  double minDistance = 1e12;
  isIncluded = true;

  std::vector<std::vector<double>> xCurrent;
  VecTool::resize(xCurrent, cell.nNodesInCell, dim);

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      for (int d = 0; d < dim; d++)
      {
        xCurrent[p][d] = cell(ic).x[p][d];
      }
    }
    dx = fabs(xCurrent[0][0] - xCurrent[1][0]);
    dy = fabs(xCurrent[0][1] - xCurrent[3][1]);
    dz = fabs(xCurrent[0][2] - xCurrent[4][2]);

    std::vector<double> N(cell.nNodesInCell, 0e0);
    ShapeFunction3D::C3D8_N(N, 0e0, 0e0, 0e0);

    for (int d = 0; d < dim; d++)
    {
      point[d] = 0e0;
      for (int p = 0; p < cell.nNodesInCell; p++)
      {
        point[d] += N[p] * node.x[cell(ic).node[p]][d];
      }
    }

    distance = 0e0;
    for (int d = 0; d < dim; d++)
    {
      diff[d] = fabs(point[d] - center[d]);
      distance += diff[d] * diff[d];
    }

    distance = sqrt(distance);
    if (distance < minDistance)
    {
      minDistance = distance;
      centerCell = ic;
      for (int d = 0; d < dim; d++)
      {
        minDiff[d] = diff[d];
      }
    }
  }

  if (minDiff[0] > dx / 2e0)
    isIncluded = false;
  if (minDiff[1] > dy / 2e0)
    isIncluded = false;
  if (minDiff[2] > dz / 2e0)
    isIncluded = false;
}

void VoxelInfo::interpolate(Node &node, Cell &cell, std::vector<std::vector<double>> &_v,
                            const int &t, const int &dim)
{
  if (isIncluded == false)
    return;

  std::vector<std::vector<double>> xCurrent;
  VecTool::resize(xCurrent, cell.nNodesInCell, dim);

  for (int p = 0; p < cell.nNodesInCell; p++)
  {
    for (int d = 0; d < dim; d++)
    {
      xCurrent[p][d] = cell(centerCell).x[p][d];
    }
  }

  double dx = fabs(xCurrent[0][0] - xCurrent[1][0]);
  double dy = dx;
  double dz = dx;
  double point[3];

  std::vector<double> N(cell.nNodesInCell, 0e0);
  ShapeFunction3D::C3D8_N(N, 0e0, 0e0, 0e0);
  for (int d = 0; d < dim; d++)
  {
    point[d] = 0e0;
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      point[d] += N[p] * xCurrent[p][d];
    }
  }

  double ss = (center[0] - point[0]);
  double tt = (center[1] - point[1]);
  double uu = (center[2] - point[2]);

  ss = ss / (dx / 2e0);
  tt = tt / (dy / 2e0);
  uu = uu / (dz / 2e0);

  if (ss < -1 - EPS || ss > 1 + EPS)
  {
    PetscPrintf(MPI_COMM_WORLD, "\ns interpolation error found.\n");
  }
  else if (tt < -1 - EPS || tt > 1 + EPS)
  {
    PetscPrintf(MPI_COMM_WORLD, "\nt interpolation error found.\n");
  }
  else if (uu < -1 - EPS || uu > 1 + EPS)
  {
    PetscPrintf(MPI_COMM_WORLD, "\nu interpolation error found.\n");
  }

  for (int p = 0; p < cell.nNodesInCell; p++)
  {
    N[p] = 0e0;
  }
  ShapeFunction3D::C3D8_N(N, ss, tt, uu);

  for (int d = 0; d < dim; d++)
  {
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      vCFD[t][d] += N[p] * _v[cell(centerCell).node[p]][d];
    }
  }
}

void VoxelInfo::average(Cell &cell, std::vector<std::vector<double>> &_v,
                        const int t, const int dim)
{
  if (cellChildren.size() == 0)
    return;
  double a = 0.1 * std::min(std::min(dx, dy), dz);

  double weightIntegral = 0e0;
  Function func3d(nNodesInCell, dim);

  for (int ic = 0; ic < cellChildren.size(); ic++)
  {
    std::vector<std::vector<double>> velCurrent;
    velCurrent.resize(nNodesInCell, std::vector<double>(dim, 0e0));

    for (int p = 0; p < nNodesInCell; p++)
    {
      for (int d = 0; d < dim; d++)
      {
        velCurrent[p][d] = _v[cell(cellChildren[ic]).node[p]][d];
        func3d.xCurrent[p][d] = cell(cellChildren[ic]).x[p][d];
      }
    }

    std::vector<double> smoothing(nNodesInCell, 1e0);
    // for(int p=0; p<nNodesInCell; p++)
    // smoothing[p] = compSmoothing(func3d, a, dim, p);

    double dxdr[3][3];
    Gauss gauss(2);

    for (int i1 = 0; i1 < 2; i1++)
    {
      for (int i2 = 0; i2 < 2; i2++)
      {
        for (int i3 = 0; i3 < 2; i3++)
        {
          ShapeFunction3D::C3D8_N(func3d.N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
          ShapeFunction3D::C3D8_dNdr(func3d.dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
          MathCommon::comp_dxdr(dxdr, func3d.dNdr, func3d.xCurrent, nNodesInCell);
          func3d.detJ = MathCommon::compDeterminant_3x3(dxdr);
          func3d.weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
          gaussIntegral(func3d, velCurrent, smoothing, weightIntegral, nNodesInCell, t, dim);
        }
      }
    }
  }
  for (int d = 0; d < dim; d++)
    vCFD[t][d] /= weightIntegral;
}

double VoxelInfo::compSmoothing(Function &func, const double a, const int dim, const int p)
{
  double sinc_x, sinc_y, sinc_z;
  double kai_x, kai_y, kai_z;

  std::vector<double> center(dim, 0e0);
  ShapeFunction3D::C3D8_N(func.N, 0e0, 0e0, 0e0);
  for (int p = 0; p < nNodesInCell; p++)
  {
    for (int d = 0; d < dim; d++)
    {
      center[d] += func.N[p] * func.xCurrent[p][d];
    }
  }

  double x0 = center[0];
  double y0 = center[1];
  double z0 = center[2];

  auto sinc = [](double x)
  {
    if (x == 0)
      return 1e0;
    return sin(PI * x) / (PI * x);
  };

  sinc_x = sinc((func.xCurrent[p][0] - x0) / dx);
  sinc_y = sinc((func.xCurrent[p][1] - y0) / dy);
  sinc_z = sinc((func.xCurrent[p][2] - z0) / dz);

  kai_x = 1 / (1 + exp(-(func.xCurrent[p][0] - (x0 - 4 * dx / 2)) / a)) - 1 / (1 + exp(-(func.xCurrent[p][0] - (x0 + 4 * dx / 2)) / a));
  kai_y = 1 / (1 + exp(-(func.xCurrent[p][1] - (y0 - 4 * dy / 2)) / a)) - 1 / (1 + exp(-(func.xCurrent[p][1] - (y0 + 4 * dy / 2)) / a));
  kai_z = 1 / (1 + exp(-(func.xCurrent[p][2] - (z0 - 4 * dz / 2)) / a)) - 1 / (1 + exp(-(func.xCurrent[p][2] - (z0 + 4 * dz / 2)) / a));

  return sinc_x * sinc_y * sinc_z * kai_x * kai_y * kai_z;
}

void VoxelInfo::gaussIntegral(Function &func, std::vector<std::vector<double>> &velCurrent, std::vector<double> &smoothing,
                              double &weightIntegral, const int nNodesInCell, const int t, const int dim)
{
  func.vol = func.detJ * func.weight;
  for (int d = 0; d < dim; d++)
  {
    for (int p = 0; p < nNodesInCell; p++)
    {
      vCFD[t][d] += func.N[p] * velCurrent[p][d] * func.vol;
    }
  }
  for (int p = 0; p < nNodesInCell; p++)
    weightIntegral += func.N[p] * smoothing[p] * func.vol;
}
