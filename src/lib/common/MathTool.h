/**
 * @file MathTool.h
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#ifndef MATHTOOL_H
#define MATHTOOL_H

#include "Array.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <functional>

using namespace std;

class MathTools2D
{
public:
  MathTools2D(const int nNodesInCell);

  double detJ, weight, vol;
  double dxdr[2][2];
  Array2D<double> xCurrent;
  Array1D<double> N;
  Array2D<double> dNdr;
  Array2D<double> dNdx;
  Array2D<double> K;

  void setZero();

  static void compInverseMatrix(double (&inv_a)[2][2], const double (&a)[2][2]);
  static double compDeterminant(const double (&a)[2][2]);

  static void comp_dxdr(double (&dxdr)[2][2], Array2D<double> &dNdr, Array2D<double> &x1, const int &nNodesInCell);
  static void comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[2][2],
                        const int &nNodesInCell);
  static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};

class MathTools3D
{
public:
  MathTools3D(const int nNodesInCell);

  double detJ, weight, vol;
  double dxdr[3][3];
  Array2D<double> xCurrent;
  Array1D<double> N;
  Array2D<double> dNdr;
  Array2D<double> dNdx;
  Array2D<double> K;

  void setZero();

  static void compInverseMatrix(double (&inv_a)[3][3], const double (&a)[3][3]);
  static double compDeterminant(const double (&a)[3][3]);

  static void comp_dxdr(double (&dxdr)[3][3], Array2D<double> &dNdr, Array2D<double> &x1, const int &nNodesInCell);
  static void comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[3][3],
                        const int &nNodesInCell);

  static double xDerivative(const std::function<double(int, int, int)> &f, int i, int j, int k, double dx);
  static double yDerivative(const std::function<double(int, int, int)> &f, int i, int j, int k, double dy);
  static double zDerivative(const std::function<double(int, int, int)> &f, int i, int j, int k, double dz);
};

#endif