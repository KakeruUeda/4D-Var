/**
 * @file MathCommon.h
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#ifndef MATHCOMMON_H
#define MATHCOMMON_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class MathCommon
{
public:
  static void compostinverseMatrix_2x2(double (&inv_a)[2][2], const double (&a)[2][2]);
  static void compostinverseMatrix_3x3(double (&inv_a)[3][3], const double (&a)[3][3]);
  static double compDeterminant_2x2(const double (&a)[2][2]);
  static double compDeterminant_3x3(const double (&a)[3][3]);
  static void comp_dxdr2D(double (&dxdr)[2][2], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
  static void comp_dNdx2D(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[2][2], const int &nNodesInCell);
  static void comp_dxdr(double (&dxdr)[3][3], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
  static void comp_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[3][3], const int &nNodesInCell);
  static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};

#endif