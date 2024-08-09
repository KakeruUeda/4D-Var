/**
 * @file PostInverseProblem.h
 * @author K.Ueda
 * @date July, 2024
 */

#include <iostream>
#include <mpi.h>
#include <map>
#include <memory>
#include <sys/stat.h>
#include "Config.h"
#include "Grid.h"
#include "MyMPI.h"
#include "VTK.h"

class PostInverseProblem
{
public:
  Application app;
  Grid grid;
  CrossSection crossSection;

  int dim;
  int nRef;
  int crossPoint;
  int flowRateVelDir;
  std::string outputDir;
  std::vector<std::vector<std::vector<double>>> velRef;
  std::vector<std::vector<std::vector<double>>> velOpt;
  void initialize(Config &conf);
  double compFlowRate(std::vector<std::vector<double>> &vel);
  double compFlowRateError(const double flowRateRefVec, const double flowRateOptVec);
};