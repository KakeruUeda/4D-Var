/**
 * @file MathTool.cpp
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#include "MathTool.h"

MathTools2D::MathTools2D(const int nNodesInCell)
{
  this->nNodesInCell = nNodesInCell;
  N.allocate(nNodesInCell);
  xCurrent.allocate(nNodesInCell, 2);
  dNdr.allocate(nNodesInCell, 2);
  dNdx.allocate(nNodesInCell, 2);
  K.allocate(nNodesInCell, nNodesInCell);
}

MathTools3D::MathTools3D(const int nNodesInCell)
{
  this->nNodesInCell = nNodesInCell;
  N.allocate(nNodesInCell);
  xCurrent.allocate(nNodesInCell, 3);
  dNdr.allocate(nNodesInCell, 3);
  dNdx.allocate(nNodesInCell, 3);
  K.allocate(nNodesInCell, nNodesInCell);
}

void MathTools2D::setZero()
{
  N.fillZero();
  xCurrent.fillZero();
  dNdr.fillZero();
  dNdx.fillZero();
  K.fillZero();
}

void MathTools3D::setZero()
{
  N.fillZero();
  xCurrent.fillZero();
  dNdr.fillZero();
  dNdx.fillZero();
  K.fillZero();
}

void MathTools2D::compInverseMatrix(double (&inv_a)[2][2], const double (&a)[2][2])
{
  double det;
  det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

  inv_a[0][0] = a[1][1] / det;
  inv_a[0][1] = -a[0][1] / det;
  inv_a[1][0] = -a[1][0] / det;
  inv_a[1][1] = a[0][0] / det;
}

void MathTools3D::compInverseMatrix(double (&inv_a)[3][3], const double (&a)[3][3])
{
  double det;
  det = compDeterminant(a);

  inv_a[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  inv_a[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
  inv_a[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
  inv_a[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
  inv_a[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
  inv_a[1][2] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
  inv_a[2][0] = a[1][0] * a[2][1] - a[1][1] * a[2][0];
  inv_a[2][1] = a[0][1] * a[2][0] - a[0][0] * a[2][1];
  inv_a[2][2] = a[0][0] * a[1][1] - a[0][1] * a[1][0];

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++)
      inv_a[i][j] = inv_a[i][j] / det;
  }
}

double MathTools2D::compDeterminant(const double (&a)[2][2])
{
  double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  return det;
}

double MathTools3D::compDeterminant(const double (&a)[3][3])
{
  double det = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2] -
               a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}

void MathTools2D::comp_dxdr(double (&dxdr)[2][2], Array2D<double> &dNdr, Array2D<double> &x1, const int numOfNodeInElm)
{
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      dxdr[i][j] = 0e0;
      for(int p = 0; p < numOfNodeInElm; p++) {
        dxdr[i][j] += dNdr(p, j) * x1(p, i);
      }
    }
  }
}

void MathTools2D::comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[2][2],
                            const int numOfNodeInElm)
{
  double drdx[2][2];
  MathTools2D::compInverseMatrix(drdx, dxdr);

  for(int p = 0; p < numOfNodeInElm; p++) {
    for(int i = 0; i < 2; i++) {
      dNdx(p, i) = 0e0;
      for(int j = 0; j < 2; j++) {
        dNdx(p, i) += dNdr(p, j) * drdx[j][i];
      }
    }
  }
}

void MathTools3D::comp_dxdr(double (&dxdr)[3][3], Array2D<double> &dNdr, Array2D<double> &x1, const int &numOfNodeInElm)
{
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      dxdr[i][j] = 0e0;
      for(int p = 0; p < numOfNodeInElm; p++) {
        dxdr[i][j] += dNdr(p, j) * x1(p, i);
      }
    }
  }
}

void MathTools3D::comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[3][3],
                            const int &numOfNodeInElm)
{
  double drdx[3][3];
  MathTools3D::compInverseMatrix(drdx, dxdr);

  for(int p = 0; p < numOfNodeInElm; p++) {
    for(int i = 0; i < 3; i++) {
      dNdx(p, i) = 0e0;
      for(int j = 0; j < 3; j++) {
        dNdx(p, i) += dNdr(p, j) * drdx[j][i];
      }
    }
  }
}

void MathTools2D::setShapesInGauss(Gauss &gauss, const int i1, const int i2)
{
  ShapeFunction2D::C2D4_N(N, gauss.point[i1], gauss.point[i2]);
  ShapeFunction2D::C2D4_dNdr(dNdr, gauss.point[i1], gauss.point[i2]);
}

void MathTools2D::setFactorsInGauss(Gauss &gauss, const int i1, const int i2)
{
  MathTools2D::comp_dxdr(dxdr, dNdr, xCurrent, nNodesInCell);
  MathTools2D::comp_dNdx(dNdx, dNdr, dxdr, nNodesInCell);
  detJ = MathTools2D::compDeterminant(dxdr);
  weight = gauss.weight[i1] * gauss.weight[i2];
  vol = detJ * weight;
}

void MathTools3D::setShapesInGauss(Gauss &gauss, const int i1, const int i2, const int i3)
{
  ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
  ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
}

void MathTools3D::setFactorsInGauss(Gauss &gauss, const int i1, const int i2, const int i3)
{
  MathTools3D::comp_dxdr(dxdr, dNdr, xCurrent, nNodesInCell);
  MathTools3D::comp_dNdx(dNdx, dNdr, dxdr, nNodesInCell);
  detJ = MathTools3D::compDeterminant(dxdr);
  weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
  vol = detJ * weight;
}

double MathTools3D::getScalarValueGP(Array1D<double> &nodeValues)
{
  double value;
  for(int p = 0; p < nNodesInCell; p++) {
    value += N(p) * nodeValues(p);
  }
  return value;
}

std::vector<double> MathTools3D::getVectorValuesGP(Array2D<double> &nodeValues)
{
  std::vector<double> values(3, 0e0);
  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < nNodesInCell; p++) {
      values[d] += N(p) * nodeValues(p, d);
    }
  }
  return values;
}