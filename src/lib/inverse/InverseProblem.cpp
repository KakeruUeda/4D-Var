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
    : main(conf), adjoint(conf), data(conf, main.grid, main.snap), app(conf.app), vel_space(conf.vel_space),
      vel_time(conf.vel_time), dim(conf.dim), outputDir(conf.outputDir), outputItr(conf.outputItr), nOMP(conf.nOMP),
      aCF(conf.aCF), bCF(conf.bCF), gCF(conf.gCF), alphaX0(conf.alphaX0), alphaX(conf.alphaX), loopMax(conf.loopMax),
      planeDir(conf.planeDir)
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
 * @brief visualize domain
 */
void InverseProblem::outputDomain()
{
  if(mpi.myId > 0) return;

  std::string vtuFile, vtiFile;

  vtuFile = outputDir + "/domain/meshPartition.vtu";
  EXPORT::exportMeshPartitionVTU(vtuFile, main.grid.node, main.grid.cell);
  vtuFile = outputDir + "/domain/phi.vtu";
  EXPORT::exportPhiVTU(vtuFile, main.grid.node, main.grid.cell);
  vtiFile = outputDir + "/domain/mask.vti";
  data.exportMaskVTI(vtiFile);
}

/**
 * @brief Run inverse routine
 */
void InverseProblem::runSimulation()
{
  outputDomain();
  compInitialOptimalVelocityField();

  // Array3D<double> v;
  // v.allocate(76, 95 * 95 * 65, 3);

  // std::string dir = "../../direct/voxelDataCreation/output/Ubend_assimilated_BC/bin";
  // for(int t = 0; t < 76; t++) {
  //   if(mpi.myId == 0) std::cout << "t = " << t << std::endl;
  //   std::string velFile = dir + "/velocityReference_" + std::to_string(t) + ".bin";
  //   v.importBIN(velFile, t);
  // }

  // for(int t = 0; t < 76; t++) {
  //   for(int i = 0; i < main.grid.nNodesGlobal; i++) {
  //     for(int d = 0; d < 3; d++) {
  //       main.vt(t, i, d) = v(t, main.grid.vecFluidUniqueNodes[i], d);
  //     }
  //   }
  // }

  // compCostFunction();

  std::ofstream cf(outputDir + "/dat/costFunction.dat");

  for(int loop = 0; loop < loopMax; loop++) {
    PetscPrintf(MPI_COMM_WORLD, "\nOpt itr. : %d\n", loop);

    if(loop == 0) compCostFunction();
    PetscPrintf(MPI_COMM_WORLD, "\ncostFunction = %e\n", costFunction.total);
    costFunction.history.push_back(costFunction.total);

    if(mpi.myId == 0) {
      cf << costFunction.term1 << " " << costFunction.term2 << " " << costFunction.term3 << " " << costFunction.term4
         << " " << costFunction.term5 << " " << costFunction.term6 << " " << costFunction.term7 << " "
         << costFunction.total;
      cf << std::endl;
    }

    bool isConverged = checkConvergence(cf, loop);
    if(isConverged) break;

    data.comp_v_err_refined();
    compFeedbackForce();

    adjoint.solveAdjoint(main, inletCB);

    compOptimalCondition();

    if(loop % outputItr == 0) {
      outputFowardSolutions(loop);
      //outputAdjointSolutions(loop);
      //outputVelocityBIN(loop);
      //outputControlVariables(loop);
      outputVelocityData(loop);
      //outputFeedbackForce(loop);
      //outputGradients(loop);
    }
    MPI_Barrier(MPI_COMM_WORLD);

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

void InverseProblem::compInletAveragedVeloicty(std::vector<std::array<double, 2>> &velArr)
{
  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double velSum = 0e0;
    double velAve = 0e0;
    int count = 0;

    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(j == 0) {
            if(fabs(data.voxel(k, j, i).v_mri(t, 1)) < 1e-8) {
              continue;
            }
            velSum += data.voxel(k, j, i).v_mri(t, 1);
            count++;
          }
        }
      }
    }

    velAve = velSum / count;

    if(mpi.myId == 0) std::cout << "t = " << t << " velAve = " << velAve << std::endl;

    std::array<double, 2> arr = {t * 0.02947812, velAve};
    velArr.push_back(arr);
  }

  MPI_Finalize();
  exit(1);
}

void InverseProblem::compInletFlowRate(std::vector<std::array<double, 2>> &velArr)
{
  double dt = main.dt * data.snap.snapInterval;

  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double flowRate = 0e0;

    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(data.voxel(k, j, i).mask < 0.5) continue;
          if(j == 0) flowRate += data.voxel(k, j, i).v_mri(t, 1) * data.dxData * data.dzData;
        }
      }
    }

    if(mpi.myId == 0) {
      std::cout << "flowRate = " << flowRate << std::endl;
    }

    std::array<double, 2> arr = {t * dt, flowRate};
    velArr.push_back(arr);
  }

  // MPI_Finalize();
  // exit(1);
}

void InverseProblem::compInletFlowRate_X(std::vector<std::array<double, 2>> &velArr)
{
  double dt = main.dt * data.snap.snapInterval;

  int sub_nx = 500;
  int sub_nz = 500;

  int nSubCell = sub_nx * sub_nz;
  int nSubNode = (sub_nx + 1) * (sub_nz + 1);

  double sub_dx = main.grid.lx / sub_nx;
  double sub_dz = main.grid.lz / sub_nz;

  std::vector<std::vector<double>> sub_x(nSubNode, std::vector<double>(2, 0e0));

  for(int k = 0; k < sub_nz + 1; k++) {
    for(int i = 0; i < sub_nx + 1; i++) {
      sub_x[i + k * (sub_nx + 1)][0] = i * sub_dx;
      sub_x[i + k * (sub_nx + 1)][1] = k * sub_dz;
    }
  }

  mt2d.nNodesInCell = 4;
  mt2d.N.allocate(4);
  mt2d.xCurrent.allocate(4, 2);
  mt2d.dNdr.allocate(4, 2);
  mt2d.dNdx.allocate(4, 2);

  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double flowRate = 0e0;
    std::vector<double> vel_y(data.nxData * data.nzData, 0e0);

    int count = 0;
    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(data.voxel(k, j, i).mask < 0.5) continue;
          if(j == 0) {
            vel_y[count] = data.voxel(k, j, i).v_mri(t, 1);
            count++;
          }
        }
      }
    }

    std::vector<double> vel_y_int(nSubNode, 0e0);

    for(int k = 0; k < sub_nz + 1; k++) {
      for(int i = 0; i < sub_nx + 1; i++) {
        int in = i + k * (sub_nx + 1);

        double px = i * sub_dx - data.dxData / 2e0;
        double pz = k * sub_dz - data.dzData / 2e0;

        int ix = px / data.dxData + EPS;
        int iz = pz / data.dzData + EPS;

        if(px < 0 || px >= main.grid.lx - data.dxData || pz < 0 || pz >= main.grid.lz - data.dzData) {
          continue;
        }

        double s = (px - (ix * data.dxData + (data.dxData / 2e0))) / (data.dxData / 2e0);
        double u = (pz - (iz * data.dzData + (data.dzData / 2e0))) / (data.dzData / 2e0);

        if(s < -1 - EPS || s > 1 + EPS || u < -1 - EPS || u > 1 + EPS) {
          if(mpi.myId == 0)
            std::cerr << "px = " << px << " pz = " << pz << " ix = " << ix << " iz = " << iz << " s = " << s
                      << " u = " << u << "\n";
          throw std::runtime_error("Interpolation error: s or u out of bounds.");
        }

        Array1D<double> N(4);
        ShapeFunction2D::C2D4_N(N, s, u);
        vel_y_int[in] = N(0) * vel_y[ix + iz * data.nxData] + N(1) * vel_y[(ix + 1) + iz * data.nxData] +
                        N(2) * vel_y[(ix + 1) + (iz + 1) * data.nxData] + N(3) * vel_y[ix + (iz + 1) * data.nxData];
      }
    }

    for(int k = 0; k < sub_nz; k++) {
      for(int i = 0; i < sub_nx; i++) {
        for(int d = 0; d < 2; d++) {
          mt2d.xCurrent(0, d) = sub_x[i + k * (sub_nx + 1)][d];
          mt2d.xCurrent(1, d) = sub_x[(i + 1) + k * (sub_nx + 1)][d];
          mt2d.xCurrent(2, d) = sub_x[(i + 1) + (k + 1) * (sub_nx + 1)][d];
          mt2d.xCurrent(3, d) = sub_x[i + (k + 1) * (sub_nx + 1)][d];
        }

        Gauss g2(2);

        for(int i2 = 0; i2 < 2; i2++) {
          for(int i1 = 0; i1 < 2; i1++) {
            mt2d.setShapesInGauss(g2, i1, i2);
            mt2d.setFactorsInGauss(g2, i1, i2);

            flowRate += mt2d.vol * (vel_y_int[i + k * (sub_nx + 1)] * mt2d.N(0) +
                                    vel_y_int[(i + 1) + k * (sub_nx + 1)] * mt2d.N(1) +
                                    vel_y_int[(i + 1) + (k + 1) * (sub_nx + 1)] * mt2d.N(2) +
                                    vel_y_int[i + (k + 1) * (sub_nx + 1)] * mt2d.N(3));
          }
        }
      }
    }

    if(mpi.myId == 0) {
      std::cout << "flowRate = " << flowRate << std::endl;
    }

    std::array<double, 2> arr = {t * dt, flowRate};
    velArr.push_back(arr);
  }
}

void InverseProblem::compInletFlowRate_top_surface(std::vector<std::array<double, 2>> &velArr)
{
  double dt = main.dt * data.snap.snapInterval;

  int sub_nx = 500;
  int sub_nz = 500;

  int nSubCell = sub_nx * sub_nz;
  int nSubNode = (sub_nx + 1) * (sub_nz + 1);

  double sub_dx = main.grid.lx / sub_nx;
  double sub_dz = main.grid.lz / sub_nz;

  std::vector<std::vector<double>> sub_x(nSubNode, std::vector<double>(2, 0e0));

  for(int k = 0; k < sub_nz + 1; k++) {
    for(int i = 0; i < sub_nx + 1; i++) {
      sub_x[i + k * (sub_nx + 1)][0] = i * sub_dx;
      sub_x[i + k * (sub_nx + 1)][1] = k * sub_dz;
    }
  }

  mt2d.nNodesInCell = 4;
  mt2d.N.allocate(4);
  mt2d.xCurrent.allocate(4, 2);
  mt2d.dNdr.allocate(4, 2);
  mt2d.dNdx.allocate(4, 2);

  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double flowRate = 0e0;
    std::vector<double> vel_y(data.nxData * data.nzData, 0e0);

    int count = 0;
    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(data.voxel(k, j, i).mask < 0.5) continue;
          if(j == data.nyData - 1) {
            vel_y[count] = data.voxel(k, j, i).v_mri(t, 1);
            count++;
          }
        }
      }
    }

    std::vector<double> vel_y_int(nSubNode, 0e0);

    for(int k = 0; k < sub_nz + 1; k++) {
      for(int i = 0; i < sub_nx + 1; i++) {
        int in = i + k * (sub_nx + 1);

        double px = i * sub_dx - data.dxData / 2e0;
        double pz = k * sub_dz - data.dzData / 2e0;

        int ix = px / data.dxData + EPS;
        int iz = pz / data.dzData + EPS;

        if(px < 0 || px >= main.grid.lx - data.dxData || pz < 0 || pz >= main.grid.lz - data.dzData) {
          continue;
        }

        double s = (px - (ix * data.dxData + (data.dxData / 2e0))) / (data.dxData / 2e0);
        double u = (pz - (iz * data.dzData + (data.dzData / 2e0))) / (data.dzData / 2e0);

        if(s < -1 - EPS || s > 1 + EPS || u < -1 - EPS || u > 1 + EPS) {
          if(mpi.myId == 0)
            std::cerr << "px = " << px << " pz = " << pz << " ix = " << ix << " iz = " << iz << " s = " << s
                      << " u = " << u << "\n";
          throw std::runtime_error("Interpolation error: s or u out of bounds.");
        }

        Array1D<double> N(4);
        ShapeFunction2D::C2D4_N(N, s, u);
        vel_y_int[in] = N(0) * vel_y[ix + iz * data.nxData] + N(1) * vel_y[(ix + 1) + iz * data.nxData] +
                        N(2) * vel_y[(ix + 1) + (iz + 1) * data.nxData] + N(3) * vel_y[ix + (iz + 1) * data.nxData];
      }
    }

    for(int k = 0; k < sub_nz; k++) {
      for(int i = 0; i < sub_nx; i++) {
        for(int d = 0; d < 2; d++) {
          mt2d.xCurrent(0, d) = sub_x[i + k * (sub_nx + 1)][d];
          mt2d.xCurrent(1, d) = sub_x[(i + 1) + k * (sub_nx + 1)][d];
          mt2d.xCurrent(2, d) = sub_x[(i + 1) + (k + 1) * (sub_nx + 1)][d];
          mt2d.xCurrent(3, d) = sub_x[i + (k + 1) * (sub_nx + 1)][d];
        }

        Gauss g2(2);

        for(int i2 = 0; i2 < 2; i2++) {
          for(int i1 = 0; i1 < 2; i1++) {
            mt2d.setShapesInGauss(g2, i1, i2);
            mt2d.setFactorsInGauss(g2, i1, i2);

            flowRate += mt2d.vol * (vel_y_int[i + k * (sub_nx + 1)] * mt2d.N(0) +
                                    vel_y_int[(i + 1) + k * (sub_nx + 1)] * mt2d.N(1) +
                                    vel_y_int[(i + 1) + (k + 1) * (sub_nx + 1)] * mt2d.N(2) +
                                    vel_y_int[i + (k + 1) * (sub_nx + 1)] * mt2d.N(3));
          }
        }
      }
    }

    // if(flowRate < -6e-5) flowRate = -6e-5;

    if(mpi.myId == 0) {
      std::cout << "flowRate = " << flowRate << std::endl;
    }

    std::array<double, 2> arr = {t * dt, flowRate};
    velArr.push_back(arr);
  }
}

void InverseProblem::compInletFlowRate_left_surface(std::vector<std::array<double, 2>> &velArr)
{
  double dt = main.dt * data.snap.snapInterval;

  int sub_ny = 500;
  int sub_nz = 500;

  int nSubCell = sub_ny * sub_nz;
  int nSubNode = (sub_ny + 1) * (sub_nz + 1);

  double sub_dy = main.grid.ly / sub_ny;
  double sub_dz = main.grid.lz / sub_nz;

  std::vector<std::vector<double>> sub_x(nSubNode, std::vector<double>(2, 0e0));

  for(int k = 0; k < sub_nz + 1; k++) {
    for(int j = 0; j < sub_ny + 1; j++) {
      sub_x[j + k * (sub_ny + 1)][0] = j * sub_dy;
      sub_x[j + k * (sub_ny + 1)][1] = k * sub_dz;
    }
  }

  mt2d.nNodesInCell = 4;
  mt2d.N.allocate(4);
  mt2d.xCurrent.allocate(4, 2);
  mt2d.dNdr.allocate(4, 2);
  mt2d.dNdx.allocate(4, 2);

  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double flowRate = 0e0;
    std::vector<double> vel_x(data.nyData * data.nzData, 0e0);

    int count = 0;
    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(data.voxel(k, j, i).mask < 0.5) continue;
          if(i == 0) {
            vel_x[count] = data.voxel(k, j, i).v_mri(t, 0);
            count++;
          }
        }
      }
    }

    std::vector<double> vel_x_int(nSubNode, 0e0);

    for(int k = 0; k < sub_nz + 1; k++) {
      for(int j = 0; j < sub_ny + 1; j++) {
        int in = j + k * (sub_ny + 1);

        double py = j * sub_dy - data.dyData / 2e0;
        double pz = k * sub_dz - data.dzData / 2e0;

        int iy = py / data.dyData + EPS;
        int iz = pz / data.dzData + EPS;

        if(py < 0 || py >= main.grid.ly - data.dyData || pz < 0 || pz >= main.grid.lz - data.dzData) {
          continue;
        }

        double s = (py - (iy * data.dyData + (data.dyData / 2e0))) / (data.dyData / 2e0);
        double u = (pz - (iz * data.dzData + (data.dzData / 2e0))) / (data.dzData / 2e0);

        if(s < -1 - EPS || s > 1 + EPS || u < -1 - EPS || u > 1 + EPS) {
          if(mpi.myId == 0)
            std::cerr << "py = " << py << " pz = " << pz << " iy = " << iy << " iz = " << iz << " s = " << s
                      << " u = " << u << "\n";
          throw std::runtime_error("Interpolation error: s or u out of bounds.");
        }

        Array1D<double> N(4);
        ShapeFunction2D::C2D4_N(N, s, u);
        vel_x_int[in] = N(0) * vel_x[iy + iz * data.nyData] + N(1) * vel_x[(iy + 1) + iz * data.nyData] +
                        N(2) * vel_x[(iy + 1) + (iz + 1) * data.nyData] + N(3) * vel_x[iy + (iz + 1) * data.nyData];
      }
    }

    for(int k = 0; k < sub_nz; k++) {
      for(int j = 0; j < sub_ny; j++) {
        for(int d = 0; d < 2; d++) {
          mt2d.xCurrent(0, d) = sub_x[j + k * (sub_ny + 1)][d];
          mt2d.xCurrent(1, d) = sub_x[(j + 1) + k * (sub_ny + 1)][d];
          mt2d.xCurrent(2, d) = sub_x[(j + 1) + (k + 1) * (sub_ny + 1)][d];
          mt2d.xCurrent(3, d) = sub_x[j + (k + 1) * (sub_ny + 1)][d];
        }

        Gauss g2(2);

        for(int i2 = 0; i2 < 2; i2++) {
          for(int i1 = 0; i1 < 2; i1++) {
            mt2d.setShapesInGauss(g2, i1, i2);
            mt2d.setFactorsInGauss(g2, i1, i2);

            flowRate += mt2d.vol * (vel_x_int[j + k * (sub_ny + 1)] * mt2d.N(0) +
                                    vel_x_int[(j + 1) + k * (sub_ny + 1)] * mt2d.N(1) +
                                    vel_x_int[(j + 1) + (k + 1) * (sub_ny + 1)] * mt2d.N(2) +
                                    vel_x_int[j + (k + 1) * (sub_ny + 1)] * mt2d.N(3));
          }
        }
      }
    }

    if(mpi.myId == 0) {
      std::cout << "flowRate = " << flowRate << std::endl;
    }

    std::array<double, 2> arr = {t * dt, flowRate};
    velArr.push_back(arr);
  }
}

void InverseProblem::compInletFlowRate_tmp(const int n, std::vector<std::array<double, 2>> &velArr)
{
  double dt = main.dt * data.snap.snapInterval;

  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double flowRate = 0e0;

    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(j == 0) {
            if(fabs(data.voxel(k, j, i).v_mri(t, 1)) > 0.5) continue;
            flowRate += data.voxel(k, j, i).v_mri(t, 1) * data.dxData * data.dzData;
          }
        }
      }
    }

    std::array<double, 2> arr;

    if(n == 0) {
      arr = {t * dt, flowRate};
      if(mpi.myId == 0) std::cout << 0.6 * flowRate << std::endl;
    } else if(n == 1) {
      arr = {t * dt, flowRate};
      if(mpi.myId == 0) std::cout << flowRate << std::endl;
    } else if(n == 2) {
      arr = {t * dt, flowRate};
    }

    velArr.push_back(arr);
  }
}

void InverseProblem::compInletMaxVelocity(std::vector<std::array<double, 2>> &velArr)
{
  double dt = main.dt * data.snap.snapInterval;

  for(int t = 0; t < data.snap.nSnapShot; t++) {
    double max = 0e0;

    for(int k = 0; k < data.nzData; k++) {
      for(int j = 0; j < data.nyData; j++) {
        for(int i = 0; i < data.nxData; i++) {
          if(j == 0) {
            if((i >= 7 && i <= 18) && (k >= 8 && k <= 19)) {
              if(fabs(data.voxel(k, j, i).v_mri(t, 1)) > 0.5) continue;
              double max_tmp = data.voxel(k, j, i).v_mri(t, 1);
              if(max_tmp > max) {  // > 0 only
                max = max_tmp;
              }
            }
          }
        }
      }
    }

    if(mpi.myId == 0) std::cout << max << std::endl;

    std::array<double, 2> arr = {t * dt, max};
    velArr.push_back(arr);
  }
}

void InverseProblem::spline_interpolate_v_mri_slice()
{
  data.u_interpolated_slice.resize(data.n_cfd_step);
  data.v_interpolated_slice.resize(data.n_cfd_step);
  data.w_interpolated_slice.resize(data.n_cfd_step);

  for(int t = 0; t < data.n_cfd_step; t++) {
    data.u_interpolated_slice[t].resize(data.u_mri_slice[0].size());
    data.v_interpolated_slice[t].resize(data.u_mri_slice[0].size());
    data.w_interpolated_slice[t].resize(data.u_mri_slice[0].size());
  }

  for(int i = 0; i < data.u_mri_slice[0].size(); i++) {
    std::vector<double> x, y1, y2, y3;

    for(int t = 0; t < data.n_mri_step; t++) {
      double p = t * data.dt_mri;
      x.push_back(p);
      y1.push_back(data.u_mri_slice[t][i]);
      y2.push_back(data.v_mri_slice[t][i]);
      y3.push_back(data.w_mri_slice[t][i]);
    }

    vector<Spline::Coefficients> cf_x = Spline::compCoefficients(x, y1);
    vector<Spline::Coefficients> cf_y = Spline::compCoefficients(x, y2);
    vector<Spline::Coefficients> cf_z = Spline::compCoefficients(x, y3);

    for(int t = 0; t < data.n_cfd_step; t++) {
      double p = t * data.dt_cfd;
      data.u_interpolated_slice[t][i] = Spline::evaluate(cf_x, p);
      data.v_interpolated_slice[t][i] = Spline::evaluate(cf_y, p);
      data.w_interpolated_slice[t][i] = Spline::evaluate(cf_z, p);
    }
  }
}

/**
 * @brief Guess initial condition.
 *        Getting initial velocity field for inverse problem
 *        using poiseuille inlet condition.
 */
void InverseProblem::compInitialOptimalVelocityField()
{
  // ----- Based On Polynomial BCs ----------------
  // spline_interpolate_v_mri_slice();
  // main.solve_NavierStokes_polynmial_BCs(data, inletCB, 3 * main.timeMax);

  // std::vector<std::map<int, std::vector<double>>> vectorVelocitySet;
  // main.solve_NavierStokes_polynmial_BCs(data, inletCB, vectorVelocitySet);
  // main.solve_NavierStokes_polynmial_BCs(data, inletCB, 0, vectorVelocitySet);
  // main.solve_NavierStokes_polynmial_BCs(data, inletCB, 1, vectorVelocitySet);
  // main.solve_NavierStokes_polynmial_BCs(data, inletCB, 2, vectorVelocitySet);

  // for(int t = 0; t < main.timeMax; t++) {
  //   for(auto &[idx, vec] : vectorVelocitySet[t]) {
  //     for(int d = 0; d < 3; d++) {
  //       X(t, idx, d) = vec[d];
  //     }
  //   }
  // }

  // for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
  //   for(int d = 0; d < 3; d++) {
  //     X0(in, d) = main.vt(0, in, d);
  //   }
  // }
  // ---------------------------------------------

  // ----- Based On Flow Rate ----------------
  std::vector<std::array<double, 2>> velArr;
  compInletFlowRate(velArr);

  // compInletFlowRate_left_surface(velArr);
  compInletFlowRate_top_surface(velArr);

  // main.solveNavierStokes(3 * main.timeMax, velArr);

  std::vector<std::map<int, std::vector<double>>> vectorVelocitySet;
  main.solveNavierStokes(velArr, vectorVelocitySet);

  for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < 3; d++) {
      main.v0(in, d) = main.vt(0, in, d);
    }
  }

  for(int t = 0; t < main.timeMax; t++) {
    for(auto &[idx, vec] : vectorVelocitySet[t]) {
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
  // ----------------------------------------
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

  for(int t = 0; t < data.n_mri_step; t++) {
    for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
      for(int d = 0; d < dim; d++) {
        data.voxel(iv).v_cfd(t, d) = 0e0;
        data.voxel(iv).v_err(t, d) = 0e0;
      }
    }
  }

  for(int t = 0; t < data.n_cfd_step; t++) {
    for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
      for(int d = 0; d < dim; d++) {
        data.voxel(iv).v_cfd_refined(t, d) = 0e0;
        data.voxel(iv).v_err_refined(t, d) = 0e0;
      }
    }
  }

  PetscPrintf(MPI_COMM_WORLD, "\nComputing Voxel Velocity ...\n");

  // ----- comp numerical velocity (space) ------
  switch(data.vel_space) {
  case VoxelVelocity::AVERAGE:
    for(int t = 0; t < data.n_cfd_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 1e-12) continue;
        if(data.voxel(iv).subId != mpi.myId) continue;
        data.average_space(main.vt, iv, t);
      }
    }
    break;

  case VoxelVelocity::WEIGHTED_AVERAGE:
    for(int t = 0; t < data.n_cfd_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 0.5) continue;
        if(data.voxel(iv).subId != mpi.myId) continue;
        data.weighted_average_space(main.vt, iv, t);
      }
    }
    break;

  case VoxelVelocity::INTERPOLATION:
    for(int t = 0; t < data.n_cfd_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 1e-12) continue;
        if(data.voxel(iv).subId != mpi.myId) continue;
        data.linear_interpolate_space(main.vt, iv, t);
      }
    }
    break;

  default:
    PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
    exit(1);
  }
  // -----------------------

  // ---- MPI communication ----
  data.gather_voxel_info();

  // ----- comp numerical velocity (time) ------
  switch(data.vel_time) {
  case VoxelVelocity::AVERAGE:
    for(int t = 0; t < data.n_mri_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        //data.average_time(iv, t);
      }
    }
    break;

  case VoxelVelocity::WEIGHTED_AVERAGE:
    for(int t = 0; t < data.n_mri_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 0.5) continue;
        data.weighted_average_time(iv, t);
      }
    }
    break;

  case VoxelVelocity::INTERPOLATION:
    for(int t = 0; t < data.n_mri_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 1e-12) continue;
        data.spline_interpolate_time(iv, t);
      }
    }
    break;

  default:
    PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
    exit(1);
  }
  // -----------------------

  // term 1
  for(int t = 0; t < data.n_mri_step; t++) {
    for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
      double dev = 0e0;
      double aCF_tmp = aCF;
      double dv_mri = data.dxData * data.dyData * data.dzData;

      for(int d = 0; d < dim; d++) {
        data.voxel(iv).v_err(t, d) = data.voxel(iv).v_cfd(t, d) - data.voxel(iv).v_mri(t, d);
        dev += data.voxel(iv).v_err(t, d) * data.voxel(iv).v_err(t, d);
      }
      if(data.voxel(iv).mask < 0.5) aCF_tmp = 0e0;
      costFunction.term1 += 5e-1 * aCF_tmp * dev * dv_mri * data.dt_mri;
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
        int n = inletCB.CBNodeMapInCell[ic][p];
        for(int d = 0; d < 2; d++) {
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
          //RegTerm2_inGaussIntegral(value2, nc, ic, t);
          RegTerm3_inGaussIntegral(value3, nc, ic, t);
          RegTerm4_inGaussIntegral(value4, nc, ic, t);
          //RegTerm5_inGaussIntegral(value5, nc, ic, t);
        }
      }
      //costFunction.term2 += 5e-1 * bCF * value2 * data.dt_cfd;
      costFunction.term3 += 5e-1 * bCF * value3 * data.dt_cfd;
      costFunction.term4 += 5e-1 * bCF * value4 * data.dt_cfd;
      //costFunction.term5 += 5e-1 * bCF * value5 * data.dt_cfd;
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
          //RegTerm6_inGaussIntegral(value6, ic);
          RegTerm7_inGaussIntegral(value7, ic);
        }
      }
    }
    //costFunction.term6 += 5e-1 * gCF * value6;
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
  mt3d.nNodesInCell = main.grid.cell.nNodesInCell;
  mt3d.N.allocate(main.grid.cell.nNodesInCell);
  mt3d.dNdr.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.xCurrent.allocate(main.grid.cell.nNodesInCell, 3);

  for(int t = 0; t < main.timeMax; t++) {
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
    if(ix < 0 || iy < 0 || iz < 0 || ix >= data.nxData || iy >= data.nyData || iz >= data.nzData) {
      return 0e0;
    }
    if(data.voxel(iz, iy, ix).mask < 0.5) {
      return 0e0;
    }
    return data.voxel(iz, iy, ix).v_err_refined(t, d);
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
 * @brief Decide step 
 * length for X0.
 */
void InverseProblem::armijoCriteriaX0(const double fk)
{
  if(alphaX0 < 1e-7) return;

  const double c1 = 1e-4;
  Array2D<double> X0_tmp;

  X0_tmp.allocate(main.grid.node.nNodesGlobal, 3);
  X0_tmp.fillZero();

  while(true) {
    if(alphaX0 < 1e-7) break;

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
  if(alphaX < 1e-7) return;

  const double c1 = 1e-4;
  Array3D<double> X_tmp;

  X_tmp.allocate(main.timeMax, main.grid.node.nNodesGlobal, 3);
  X_tmp.fillZero();

  while(true) {
    if(alphaX < 1e-7) break;

    X_tmp = X + alphaX * (-gradX);

    main.solveNavierStokes(X0, X_tmp);
    compCostFunction();

    double lk = costFunction.total;

    double tmp = 0e0;

    for(int t = 0; t < adjoint.timeMax; t++) {
      //if(t % main.snap.snapInterval == 0) {
      for(int in = 0; in < main.grid.nNodesGlobal; in++) {
        for(int d = 0; d < 3; d++) {
          tmp += -(gradX(t, in, d) * gradX(t, in, d));
        }
      }
      //}
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