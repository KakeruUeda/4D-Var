/**
 * @file    FEM.h
 * @author  K.Ueda
 * @date    August, 2024
 */

#include "FEM.h"

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

double FEM::comp_pulse(const int t)
{
  double timePhase = (t - pulseBeginItr) * dt;
  return pulse = 0.25 * sin((2e0 * PI / T) * timePhase) + 1.0;
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