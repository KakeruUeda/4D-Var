/**
 * @file    FEM.h
 * @author  K.Ueda
 * @date    August, 2024
 */

#include "FEM.h"

/**********************
 * @brief Constructor.
 */
FEM::FEM(Config &conf)
    : dim(conf.dim), nOMP(conf.nOMP), rho(conf.rho), mu(conf.mu), dt(conf.dt), timeMax(conf.timeMax),
      pulsatileFlow(conf.pulsatileFlow), pulseBeginItr(conf.pulseBeginItr), T(conf.T), alpha(conf.alpha),
      resistance(conf.resistance)
{
}

/********************************
 * @brief Compute element length.
 */
double FEM::comp_he(Array2D<double> &x)
{
  return fabs(x(1, 0) - x(0, 0));
}

/**********************************
 * @brief Compute Darcy resistance.
 */
double FEM::comp_f(const double phi)
{
  double coeff = resistance * alpha;
  double value = coeff * (1e0 - phi) / (alpha + phi);
  return coeff * value;
}

/*****************************************
 * @brief Compute stabilization parameter.
 */
double FEM::comp_tau(std::vector<double> &vel, const double he)
{
  double tau = 0e0;
  double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

  double term1 = (2e0 / dt) * (2e0 / dt);
  double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
  double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));

  return tau = pow(term1 + term2 + term3, -5e-1);
}

/***********************
 * @brief Compute pulse.
 */
double FEM::comp_tau(double vel[3], const double he)
{
  double tau = 0e0;
  double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

  double term1 = (2e0 / dt) * (2e0 / dt);
  double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
  double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));

  return tau = pow(term1 + term2 + term3, -5e-1);
}

double FEM::comp_pulse(const int t)
{
  double timePhase = (t - pulseBeginItr) * dt;
  return pulse = 0.25 * sin((2e0 * PI / T) * timePhase) + 1.0;
}

/***********************************************
 * @brief Compute arbitary value on gauss point.
 */
void FEM::compValueOnGaussPoint(Grid &grid, MathTools3D &tools, double (&arr)[3], Array2D<double> &value, const int ic)
{
  for(int d = 0; d < 3; d++) {
    arr[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      arr[d] += tools.N(p) * value(n, d);
    }
  }
}

/****************************************************
 * @brief Compute arbitary derivative on gauss point.
 */
void FEM::compDerivOnGaussPoint(Grid &grid, MathTools3D &tools, double (&arr)[3][3], Array2D<double> &value,
                                const int ic)
{
  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      arr[d1][d2] = 0e0;
      for(int p = 0; p < grid.cell.nNodesInCell; p++) {
        int n = grid.cell(ic).node[p];
        arr[d1][d2] += tools.dNdx(p, d2) * value(n, d1);
      }
    }
  }
}

/*****************************************
 * @brief Compute velocity on gauss point.
 */
void FEM::compVelOnGaussPoint(Grid &grid, MathTools3D &tools, double (&arr)[3], const int ic)
{
  for(int d = 0; d < 3; d++) {
    arr[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      arr[d] += tools.N(p) * v(n, d);
    }
  }
}

/******************************************
 * @brief Compute advection on gauss point.
 */
void FEM::compAdVelOnGaussPoint(Grid &grid, MathTools3D &tools, double (&arr)[3], const int ic)
{
  for(int d = 0; d < 3; d++) {
    arr[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      arr[d] += tools.N(p) * (1.5 * v(n, d) - 0.5 * vPrev(n, d));
    }
  }
}

/****************************************************
 * @brief Compute velocity derivative on gauss point.
 */
void FEM::compVelDerivOnGaussPoint(Grid &grid, MathTools3D &tools, double (&arr)[3][3], const int ic)
{
  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      arr[d1][d2] = 0e0;
      for(int p = 0; p < grid.cell.nNodesInCell; p++) {
        int n = grid.cell(ic).node[p];
        arr[d1][d2] += tools.dNdx(p, d2) * v(n, d1);
      }
    }
  }
}

/**************************
 * @brief Update row index.
 */
void FEM::updateRowIndex(Grid &grid, const int ii, const int ic)
{
  IU = grid.cell(ic).dofStart[ii];
  IV = IU + 1;
  IW = IU + 2;
  IP = IU + 3;
  ILU = IU + 4;
  ILV = IU + 5;
  ILW = IU + 6;
}

/*****************************
 * @brief Update column index.
 */
void FEM::updateColumnIndex(Grid &grid, const int jj, const int ic)
{
  JU = grid.cell(ic).dofStart[jj];
  JV = JU + 1;
  JW = JU + 2;
  JP = JU + 3;
  JLU = JU + 4;
  JLV = JU + 5;
  JLW = JU + 6;
}

/*****************************
 * @brief Update 2D row index.
 */
void FEM::updateRowIndexPlane(Grid &grid, const int ii, const int ic)
{
  IU = grid.cell(ic).dofStartPlane[ii];
  IV = IU + 1;
  IW = IU + 2;
  IP = IU + 3;
  ILU = IU + 4;
  ILV = IU + 5;
  ILW = IU + 6;
}

/********************************
 * @brief Update 2D column index.
 */
void FEM::updateColumnIndexPlane(Grid &grid, const int jj, const int ic)
{
  JU = grid.cell(ic).dofStartPlane[jj];
  JV = JU + 1;
  JW = JU + 2;
  JP = JU + 3;
  JLU = JU + 4;
  JLV = JU + 5;
  JLW = JU + 6;
}