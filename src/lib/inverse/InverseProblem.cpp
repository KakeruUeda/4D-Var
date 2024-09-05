/**
 * @file InverseProblem.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/**
 * @brief construct inverse object from config param
 */
InverseProblem::InverseProblem(Config &conf)
    : main(conf), adjoint(conf), data(conf, main.grid, main.snap), app(conf.app), vvox(conf.vvox), dim(conf.dim),
      outputDir(conf.outputDir), outputItr(conf.outputItr), nOMP(conf.nOMP), aCF(conf.aCF), bCF(conf.bCF),
      gCF(conf.gCF), alphaX0(conf.alphaX0), alphaX(conf.alphaX), loopMax(conf.loopMax), planeDir(conf.planeDir)
{
  std::string dir;
  std::string output = "output";

  mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  outputDir = "output/" + outputDir;
  mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/domain";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/data";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/dat";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/main";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/adjoint";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/other";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/bin";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/optimized";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

  main.outputDir = outputDir;
  adjoint.outputDir = outputDir;
}

/**
 * @brief Run inverse routine
 */
void InverseProblem::runSimulation()
{
  main.outputDomain();
  compInitialOptimalVelocityField();
  std::ofstream cf(outputDir + "/dat/costFunction.dat");

  for(int loop = 0; loop < loopMax; loop++) {
    PetscPrintf(MPI_COMM_WORLD, "\nOpt itr. : %d\n", loop);

    compCostFunction();
    PetscPrintf(MPI_COMM_WORLD, "\ncostFunction = %e\n", costFunction.total);

    costFunction.history.push_back(costFunction.total);

    if(mpi.myId == 0) {
      cf << costFunction.term1 << " " << costFunction.term2 << " " 
         << costFunction.term3 << " " << costFunction.term4 << " " 
         << costFunction.term5 << " " << costFunction.term6 << " " 
         << costFunction.term7 << " " << costFunction.total;
      cf << std::endl;
    }

    bool isConverged = checkConvergence(cf, loop);
    if(isConverged) break;

    compFeedbackForce();
    compTimeInterpolatedFeedbackForce();

    adjoint.solveAdjoint(main, inletCB);

    if(loop % outputItr == 0) {
      outputFowardSolutions(loop);
      outputAdjointSolutions(loop);
      outputVelocityBIN(loop);
      outputControlVariables(loop);
      outputVelocityData(loop);
      outputFeedbackForce(loop);
      outputGradients(loop);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    compOptimalCondition();

    armijoCriteriaX0(costFunction.total);
    armijoCriteriaX(costFunction.total);

    PetscPrintf(MPI_COMM_WORLD, "\n(alphaX0, alphaX) = (%f, %f)\n", alphaX0, alphaX);
  }
  cf.close();

  outputOptimizedVariables();
}

/**
 * @brief Check cnvergence of the cost function.
 */
bool InverseProblem::checkConvergence(std::ofstream &cf, const int loop)
{
  if(loop >= 10) {
    double diff = costFunction.history[loop] - costFunction.history[loop - 10];
    diff = fabs(diff);

    if(diff < 1e-3) {
      PetscPrintf(MPI_COMM_WORLD, "Converged. OptItr = %d", loop);
      cf.close();
      return true;
    }
  }
  return false;
}

/**
 * @brief Guess initial condition.
 *        Getting initial velocity field for inverse problem
 *        using poiseuille inlet condition.
 */
void InverseProblem::compInitialOptimalVelocityField()
{
  main.solveNaveirStokes(main.timeMax * 3);

  for(int t = 0; t < main.timeMax; t++) {
    for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
      for(int d = 0; d < 3; d++) {
        main.vt(t, in, d) = main.v0(in, d);
      }
    }
  }

  // take snapshot
  int snapCount = 0;
  for(int t = 0; t < main.timeMax; t++) {
    if((t - main.snap.snapTimeBeginItr) % main.snap.snapInterval == 0) {
      main.snap.takeSnapShot(main.vt, main.grid.node.nNodesGlobal, snapCount, t);
      snapCount++;
    }
  }

  for(int t = 0; t < main.timeMax; t++) {
    for(auto &[idx, vec] : main.dirichlet.velocitySet) {
      for(int d = 0; d < 3; d++) {
        X(t, idx, d) = vec[d];
      }
    }
  }

  for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < main.dim; d++) {
      X0(in, d) = main.v0(in, d);
    }
  }

  /*
  for(int icb = 0; icb < inletCB.CBNodeMap.size(); icb++){
    int in = inletCB.CBNodeMap[icb];
    for(int d = 0; d < 3; d++){
      X0(in, d) = 0e0;
    }
  }
  */
}

/**
 * @brief Compute discrepancies over domain and time
 */
void InverseProblem::compCostFunction()
{
  costFunction.term1 = 0e0;
  costFunction.term2 = 0e0;
  costFunction.term3 = 0e0;
  costFunction.term4 = 0e0;
  costFunction.term5 = 0e0;
  costFunction.term6 = 0e0;
  costFunction.term7 = 0e0;

  for(int t = 0; t < main.snap.nSnapShot; t++) {
    for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
      for(int d = 0; d < dim; d++) {
        data.voxel(iv).v_cfd(t, d) = 0e0;
        data.voxel(iv).v_err(t, d) = 0e0;
      }
    }
  }

  switch(data.vvox) {
  case VoxelVelocity::AVERAGE:
    for(int t = 0; t < main.snap.nSnapShot; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        data.average(iv, t);
      }
    }
    break;

  case VoxelVelocity::INTERPOLATION:
    for(int t = 0; t < main.snap.nSnapShot; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        //data.interpolate(iv, t);
      }
    }
    break;

  default:
    PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
    exit(1);
  }

  // term 1
  for(int t = 0; t < main.snap.nSnapShot; t++) {
    for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
      if(data.voxel(iv).mask < 0.999) {
        continue;
      }
      double dev = 0e0;
      double volume = data.dxData * data.dyData * data.dzData;
      double deltaT = main.dt * main.snap.snapInterval;
      for(int d = 0; d < dim; d++) {
        data.voxel(iv).v_err(t, d) = data.voxel(iv).v_cfd(t, d) - data.voxel(iv).v_mri(t, d);
        dev += data.voxel(iv).v_err(t, d) * data.voxel(iv).v_err(t, d);
      }
      costFunction.term1 += 5e-1 * aCF * dev * volume * deltaT;
    }
  }

  int nc = inletCB.CBNodeMapInCell[0].size();

  Gauss gauss(2);

  mt2d.nNodesInCell = nc;
  mt2d.N.allocate(nc);
  mt2d.dNdr.allocate(nc, 2);
  mt2d.dNdx.allocate(nc, 2);
  mt2d.xCurrent.allocate(nc, 2);

  mt3d.nNodesInCell = main.grid.cell.nNodesInCell;
  mt3d.N.allocate(main.grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.xCurrent.allocate(main.grid.cell.nNodesInCell, 3);

  // term 2, 3, 4, 5
  for(int t = 0; t < main.timeMax; t++) {
    for(int ic = 0; ic < inletCB.CBCellMap.size(); ic++) {

      for(int p = 0; p < nc; p++) {
        for(int d = 0; d < dim - 1; d++) {
          int n = inletCB.CBNodeMapInCell[ic][p];
          mt2d.xCurrent(p, d) = main.grid.node.x[n][planeDir[d]];
        }
      }

      double value2 = 0e0;
      double value3 = 0e0;
      double value4 = 0e0;
      double value5 = 0e0;

      for(int i1 = 0; i1 < 2; i1++) {
        for(int i2 = 0; i2 < 2; i2++) {
          mt2d.setShapesInGauss(gauss, i1, i2);
          mt2d.setFactorsInGauss(gauss, i1, i2);
          RegTerm2_inGaussIntegral(value2, nc, ic, t);
          RegTerm3_inGaussIntegral(value3, nc, ic, t);
          RegTerm4_inGaussIntegral(value4, nc, ic, t);
          RegTerm5_inGaussIntegral(value5, nc, ic, t);
        }
      }
      costFunction.term2 += 5e-1 * bCF * value2;
      costFunction.term3 += 5e-1 * bCF * value3;
      costFunction.term4 += 5e-1 * bCF * value4;
      costFunction.term5 += 5e-1 * bCF * value5;
    }
  }

  // term 6, 7
  for(int ic = 0; ic < main.grid.cell.nCellsGlobal; ic++) {

    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      for(int d = 0; d < 3; d++) {
        mt3d.xCurrent(p, d) = main.grid.node.x[main.grid.cell(ic).node[p]][d];
      }
    }

    double value6 = 0e0;
    double value7 = 0e0;

    for(int i1 = 0; i1 < 2; i1++) {
      for(int i2 = 0; i2 < 2; i2++) {
        for(int i3 = 0; i3 < 2; i3++) {
          mt3d.setShapesInGauss(gauss, i1, i2, i3);
          mt3d.setFactorsInGauss(gauss, i1, i2, i3);
          RegTerm6_inGaussIntegral(value6, ic);
          RegTerm7_inGaussIntegral(value7, ic);
        }
      }
    }
    costFunction.term6 += 5e-1 * gCF * value6;
    costFunction.term7 += 5e-1 * gCF * value7;
  }

  costFunction.sum();
}

/**
 * @brief Compute value for regularization term2
 *        on gauss integral point.
 */
void InverseProblem::RegTerm2_inGaussIntegral(double &value, const int nc, const int ic, const int t)
{
  double u[3] = {0e0, 0e0, 0e0};

  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < nc; p++) {
      int n = inletCB.CBNodeMapInCell[ic][p];
      u[d] += mt2d.N(p) * main.vt(t, n, d);
    }
  }

  for(int d = 0; d < 3; d++) {
    value += u[d] * u[d] * mt2d.vol;
  }
}

/**
 * @brief Compute value for regularization term3
 *        on gauss integral point.
 */
void InverseProblem::RegTerm3_inGaussIntegral(double &value, const int nc, const int ic, const int t)
{
  double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      for(int p = 0; p < nc; p++) {
        int n = inletCB.CBNodeMapInCell[ic][p];
        dudx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t, n, d1);
      }
    }
  }

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      value += dudx[d1][d2] * dudx[d1][d2] * mt2d.vol;
    }
  }
}

/**
 * @brief Compute value for regularization term4
 *        on gauss integral point.
 */
void InverseProblem::RegTerm4_inGaussIntegral(double &value, const int nc, const int ic, const int t)
{
  double u[3] = {0e0, 0e0, 0e0};
  double ub[3] = {0e0, 0e0, 0e0};
  double dudt[3] = {0e0, 0e0, 0e0};

  for(int p = 0; p < nc; p++) {
    int n = inletCB.CBNodeMapInCell[ic][p];
    for(int d = 0; d < 3; d++) {
      if(t == 0) {
        u[d] += mt2d.N(p) * main.vt(t, n, d);
        ub[d] += mt2d.N(p) * main.v0(n, d);
      } else {
        u[d] += mt2d.N(p) * main.vt(t, n, d);
        ub[d] += mt2d.N(p) * main.vt(t - 1, n, d);
      }
    }
  }

  for(int d = 0; d < 3; d++) {
    dudt[d] = (u[d] - ub[d]) / main.dt;
  }
  for(int d = 0; d < 3; d++) {
    value += dudt[d] * dudt[d] * mt2d.vol;
  }
}

/**
 * @brief Compute value for regularization term3
 *        on gauss integral point.
 */
void InverseProblem::RegTerm5_inGaussIntegral(double &value, const int nc, const int ic, const int t)
{
  double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  double dubdx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  double dudxdt[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  for(int p = 0; p < nc; p++) {
    int n = inletCB.CBNodeMapInCell[ic][p];
    for(int d1 = 0; d1 < 3; d1++) {
      for(int d2 = 0; d2 < 2; d2++) {
        if(t == 0) {
          dudx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t, n, d1);
          dubdx[d1][d2] += mt2d.dNdx(p, d2) * main.v0(n, d1);
        } else {
          dudx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t, n, d1);
          dubdx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t - 1, n, d1);
        }
      }
    }
  }

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      dudxdt[d1][d2] = (dudx[d1][d2] - dubdx[d1][d2]) / main.dt;
    }
  }
  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      value += dudxdt[d1][d2] * dudxdt[d1][d2] * mt2d.vol;
    }
  }
}

/**
 * @brief Compute value for regularization term6
 *        on gauss integral point.
 */
void InverseProblem::RegTerm6_inGaussIntegral(double &value, const int ic)
{
  double u0[3] = {0e0, 0e0, 0e0};

  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int n = main.grid.cell(ic).node[p];
      u0[d] += mt3d.N(p) * main.v0(n, d);
    }
  }

  for(int d = 0; d < 3; d++) {
    value += u0[d] * u0[d] * mt3d.vol;
  }
}

/**
 * @brief Compute value for regularization term7
 *        on gauss integral point.
 */
void InverseProblem::RegTerm7_inGaussIntegral(double &value, const int ic)
{
  double du0dx[3][3] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
        int n = main.grid.cell(ic).node[p];
        du0dx[d1][d2] += mt3d.dNdx(p, d2) * main.v0(n, d1);
      }
    }
  }

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      value += du0dx[d1][d2] * du0dx[d1][d2] * mt3d.vol;
    }
  }
}

/**
 * @brief Compute value for regularization term3
 *        on gauss integral point.
 */
void InverseProblem::compFeedbackForce()
{
  for(int t = 0; t < main.snap.nSnapShot; t++) {
    for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
      for(int d = 0; d < dim; d++) {
        adjoint.feedbackForce(t, in, d) = 0e0;
      }
    }
    for(int ic = 0; ic < main.grid.cell.nCellsGlobal; ic++) {
      assembleFeedbackForce(ic, t);
    }
  }
}

/**
 * @brief Assemble RHS feedback force for adjoint system.
 * @param ic: cell index.
 * @param t: time index.
 */
void InverseProblem::assembleFeedbackForce(const int ic, const int t)
{
  Gauss g2(2);

  mt3d.nNodesInCell = main.grid.cell.nNodesInCell;
  mt3d.N.allocate(main.grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.xCurrent.allocate(main.grid.cell.nNodesInCell, 3);

  for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
    for(int d = 0; d < dim; d++) {
      mt3d.xCurrent(p, d) = main.grid.cell(ic).x[p][d];
    }
  }

  for(int i1 = 0; i1 < 2; i1++) {
    for(int i2 = 0; i2 < 2; i2++) {
      for(int i3 = 0; i3 < 2; i3++) {
        mt3d.setShapesInGauss(g2, i1, i2, i3);
        mt3d.setFactorsInGauss(g2, i1, i2, i3);

        double feedback[3] = {0e0, 0e0, 0e0};
        double point[3] = {0e0, 0e0, 0e0};

        for(int d = 0; d < 3; d++) {
          for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
            point[d] += mt3d.N(p) * mt3d.xCurrent(p, d);
          }
        }
        compInterpolatedFeeback(feedback, point, t);
        feedbackGaussIntegral(feedback, ic, t);
      }
    }
  }
}

/**
 * @brief Interpolate space-discrete feedback force onto CFD node.
 */
void InverseProblem::compInterpolatedFeeback(double (&feedback)[3], double (&point)[3], const int t)
{
  double px = point[0] - (data.dxData / 2e0);
  double py = point[1] - (data.dyData / 2e0);
  double pz = point[2] - (data.dzData / 2e0);

  int ix = (px / data.dxData) + EPS;
  int iy = (py / data.dyData) + EPS;
  int iz = (pz / data.dzData) + EPS;

  double ss = (px - (ix * data.dxData + (data.dxData / 2e0)));
  double tt = (py - (iy * data.dyData + (data.dyData / 2e0)));
  double uu = (pz - (iz * data.dzData + (data.dzData / 2e0)));

  ss = ss / (data.dxData / 2e0);
  tt = tt / (data.dyData / 2e0);
  uu = uu / (data.dzData / 2e0);

  Array1D<double> N(main.grid.cell.nNodesInCell);
  ShapeFunction3D::C3D8_N(N, ss, tt, uu);

  auto velocity = [&](int iz, int iy, int ix, int d, const int t) -> double {
    if(ix < 0 || iy < 0 || iz < 0 || ix >= main.grid.nx || iy >= main.grid.ny || iz >= main.grid.nz) {
      return 0e0;
    }
    return data.voxel(iz, iy, ix).v_err(t, d);
  };

  for(int d = 0; d < 3; d++) {
    feedback[d] = N(0) * velocity(iz, iy, ix, d, t) + N(1) * velocity(iz, iy, ix + 1, d, t) +
                  N(2) * velocity(iz, iy + 1, ix + 1, d, t) + N(3) * velocity(iz, iy + 1, ix, d, t) +
                  N(4) * velocity(iz + 1, iy, ix, d, t) + N(5) * velocity(iz + 1, iy, ix + 1, d, t) +
                  N(6) * velocity(iz + 1, iy + 1, ix + 1, d, t) + N(7) * velocity(iz + 1, iy + 1, ix, d, t);
  }
}

/**
 * @brief Compute feedback force on gauss integral point.
 */
void InverseProblem::feedbackGaussIntegral(double (&feedback)[3], const int ic, const int t)
{
  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int in = main.grid.cell(ic).node[p];
      adjoint.feedbackForce(t, in, d) += aCF * mt3d.N(p) * feedback[d] * mt3d.vol;
    }
  }
}

/**
 * @brief Interpolate the feedback force at each CFD time step.
 */
void InverseProblem::compTimeInterpolatedFeedbackForce()
{
  for(int in = 0; in < adjoint.grid.node.nNodesGlobal; in++) {
    std::vector<double> x, y1, y2, y3;

    for(int t = 0; t < main.snap.nSnapShot; t++) {
      double dp = t * main.dt * main.snap.snapInterval;
      x.push_back(dp);
      y1.push_back(adjoint.feedbackForce(t, in, 0));
      y2.push_back(adjoint.feedbackForce(t, in, 1));
      y3.push_back(adjoint.feedbackForce(t, in, 2));
    }

    vector<Spline::Coefficients> cf_x = Spline::compCoefficients(x, y1);
    vector<Spline::Coefficients> cf_y = Spline::compCoefficients(x, y2);
    vector<Spline::Coefficients> cf_z = Spline::compCoefficients(x, y3);

    for(int t = 0; t < adjoint.timeMax; t++) {
      double p = t * main.dt;
      adjoint.feedbackForceT(t, in, 0) = Spline::evaluate(cf_x, p);
      adjoint.feedbackForceT(t, in, 1) = Spline::evaluate(cf_y, p);
      adjoint.feedbackForceT(t, in, 2) = Spline::evaluate(cf_z, p);
    }
  }
}


/**
 * @brief Decide step length for X0.
 */
void InverseProblem::armijoCriteriaX0(const double fk)
{
  const double c1 = 1e-4;
  Array2D<double> X0_tmp;

  X0_tmp.allocate(main.grid.node.nNodesGlobal, 3);
  X0_tmp.fillZero();

  while(true) {
    X0_tmp = X0 + alphaX0 * (-gradX0);

    main.solveNavierStokes(X0_tmp, X);
    compCostFunction();

    double lk = costFunction.total;

    double tmp = 0e0;
    for(int in = 0; in < main.grid.nNodesGlobal; in++) {
      for(int d = 0; d < 3; d++) {
        tmp += -(gradX0(in, d) * gradX0(in, d));
      }
    }

    double l_tmp = fk + c1 * tmp * alphaX0;

    if(lk <= l_tmp) {
      main.v0 = X0_tmp;
      X0 = X0_tmp;
      break;
    } else {
      alphaX0 = alphaX0 * 5e-1;
      PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alphaX0);
    }
  }
}

/**
 * @brief Decide step length for X.
 */
void InverseProblem::armijoCriteriaX(const double fk)
{
  const double c1 = 1e-4;
  Array3D<double> X_tmp;

  X_tmp.allocate(main.timeMax, main.grid.node.nNodesGlobal, 3);
  X_tmp.fillZero();

  while(true) {
    X_tmp = X + alphaX * (-gradX);

    main.solveNavierStokes(X0, X_tmp);
    compCostFunction();

    double lk = costFunction.total;

    double tmp = 0e0;
    for(int t = 0; t < adjoint.timeMax; t++) {
      if(t % main.snap.snapInterval == 0) {
        for(int in = 0; in < main.grid.nNodesGlobal; in++) {
          for(int d = 0; d < 3; d++) {
            tmp += -(gradX(t, in, d) * gradX(t, in, d));
          }
        }
      }
    }
    double l_tmp = fk + c1 * tmp * alphaX;

    if(lk <= l_tmp) {
      X = X_tmp;
      break;
    } else {
      alphaX = alphaX * 5e-1;
      PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alphaX);
    }
  }
}