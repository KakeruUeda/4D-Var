/**
 * @file MathTool.h
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#ifndef MATHTOOL_H
#define MATHTOOL_H

#include "Array.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class MathTools2D
{
public:
  MathTools2D(const int nNodesInCell);
  MathTools2D(){}

  double detJ, weight, vol;
  double dxdr[2][2];
  int nNodesInCell;
  
  Array2D<double> xCurrent;
  Array1D<double> N;
  Array2D<double> dNdr;
  Array2D<double> dNdx;
  Array2D<double> K;

  void setZero();
  void setShapesInGauss(Gauss &gauss, const int i1, const int i2);
  void setFactorsInGauss(Gauss &gauss, const int i1, const int i2);
  static void compInverseMatrix(double (&inv_a)[2][2], const double (&a)[2][2]);
  static double compDeterminant(const double (&a)[2][2]);

  static void comp_dxdr(double (&dxdr)[2][2], Array2D<double> &dNdr, Array2D<double> &x1, const int nNodesInCell);
  static void comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[2][2],
                        const int nNodesInCell);
  static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};

class MathTools3D
{
public:
  MathTools3D(const int nNodesInCell);
  MathTools3D(){}

  double detJ, weight, vol;
  double dxdr[3][3];
  int nNodesInCell;

  Array2D<double> xCurrent;
  Array1D<double> N;
  Array2D<double> dNdr;
  Array2D<double> dNdx;
  Array2D<double> K;

  void setZero();
  void setShapesInGauss(Gauss &gauss, const int i1, const int i2, const int i3);
  void setFactorsInGauss(Gauss &gauss, const int i1, const int i2, const int i3);
  static void compInverseMatrix(double (&inv_a)[3][3], const double (&a)[3][3]);
  static double compDeterminant(const double (&a)[3][3]);

  double getScalarValueGP(Array1D<double> &nodeValues);
  std::vector<double> getVectorValuesGP(Array2D<double> &nodeValues);

  static void comp_dxdr(double (&dxdr)[3][3], Array2D<double> &dNdr, Array2D<double> &x1, const int &nNodesInCell);
  static void comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[3][3],
                        const int &nNodesInCell);
};

#endif