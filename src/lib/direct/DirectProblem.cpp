/**
 * @file DirectProblem.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/**
 * @brief Constructer.
 */
DirectProblem::DirectProblem(Config &conf)
    : FEM(conf), app(conf.app), dim(conf.dim), outputDir(conf.outputDir),
      nOMP(conf.nOMP), grid(conf), snap(conf)
{
  if(app == Application::USNS) {
    std::string dir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/domain";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/data";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/solution";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/input";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  }
}

/**
 * @brief Simulate Unsteady Navier Stokes Equation.
 */
void DirectProblem::runSimulation()
{
  outputDomain();
  solveNavierStokes();
}

/**
 * @brief Visualize partitioned domain and phi.
 */
void DirectProblem::outputDomain()
{
  if(mpi.myId > 0) {
    return;
  }
  std::string vtuFile;
  std::string vtiFile;

  vtuFile = outputDir + "/domain/meshPartition.vtu";
  EXPORT::exportMeshPartitionVTU(vtuFile, grid.node, grid.cell);
  vtuFile = outputDir + "/domain/phi.vtu";
  EXPORT::exportPhiVTU(vtuFile, grid.node, grid.cell);
}
