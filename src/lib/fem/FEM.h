/**
 * @file    FEM.h
 * @author  K.Ueda
 * @date    August, 2024
 */

#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include "Array.h"
#include "Grid.h"

// FEM Interface
class FEM
{
public:
  FEM(Config &conf);

  int dim, nOMP;

  // Index
  int IU, IV, IW, IP;
  int ILU, ILV, ILW;
  int JU, JV, JW, JP;
  int JLU, JLV, JLW;

  // Time parameter
  double dt;
  int timeMax;
  int pulsatileFlow;
  int pulseBeginItr;
  double T;
  double pulse;

  // Stabilization parameter
  double tau;


  // Pysical parameter
  double Re, rho, mu, nu;

  // Darcy parameter
  double alpha, resistance;

  // Foward variable
  Array2D<double> v0;
  Array2D<double> v;
  Array2D<double> vPrev;
  Array1D<double> p;
  Array3D<double> vt;
  Array2D<double> pt;

  // Adjoint variable
  Array2D<double> w;
  Array2D<double> wPrev;
  Array2D<double> q;
  Array2D<double> qPrev;
  Array2D<double> l;
  Array3D<double> wt;
  Array3D<double> qt;
  Array3D<double> lt;

  // for vti visualization
  Array2D<double> vvti;
  Array1D<double> pvti;
  Array2D<double> wvti;
  Array1D<double> qvti;
  Array2D<double> lvti;

  double vgp[3];
  double advgp[3];
  double dvgpdx[3][3];

  double vk[3], vk1[3], vk2[3];
  double dvkdx[3][3], dvk1dx[3][3], dvk2dx[3][3];
  double advk1[3], advk2[3], advk3[3];
  double dpkdx[3], dpk1dx[3], dpk2dx[3];
  double wk[3], wk1[3], wk2[3];
  double dwkdx[3][3], dwk1dx[3][3], dwk2dx[3][3];
  double dqkdx[3], dqk1dx[3], dqk2dx[3];

  double comp_he(Array2D<double> &x);
  double comp_f(const double phi);
  double comp_tau(std::vector<double> &vel, const double he);
  double comp_pulse(const int t);

  void updateRowIndex(Grid &grid, const int ii, const int ic);
  void updateColumnIndex(Grid &grid, const int jj, const int ic);
  void updateRowIndexPlane(Grid &grid, const int ii, const int ic);
  void updateColumnIndexPlane(Grid &grid, const int jj, const int ic);
  
};


#endif
