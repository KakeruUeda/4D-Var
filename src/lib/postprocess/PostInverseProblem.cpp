#include "PostInverseProblem.h"

void PostInverseProblem::initialize(Config &conf)
{
  app = conf.app;
  dim = conf.dim;
  nRef = conf.nRef;
  flowRateVelDir = conf.flowRateVelDir;
  crossPoint = conf.crossPoint;
  velRef = conf.velRef;
  velOpt = conf.velOpt;
  outputDir = conf.outputDir;
  crossSection = conf.crossSection;
  grid.gridType = conf.gridType;

  if (grid.gridType == GridType::STRUCTURED)
  {
    grid.nx = conf.nx;
    grid.ny = conf.ny;
    grid.nz = conf.nz;
    grid.lx = conf.lx;
    grid.ly = conf.ly;
    grid.lz = conf.lz;
    grid.dx = conf.dx;
    grid.dy = conf.dy;
    grid.dz = conf.dz;

    grid.node.nNodesGlobal = (grid.nx + 1) * (grid.ny + 1) * (grid.nz + 1);
    grid.cell.nCellsGlobal = grid.nx * grid.ny * grid.nz;
    grid.cell.nNodesInCell = conf.nNodesInCell;

    if (grid.node.nNodesGlobal != velRef[0].size() || grid.node.nNodesGlobal != velOpt[0].size())
    {
      std::cout << "Num of structured grid nodes doesn't match with the num of input velocity" << std::endl;
      exit(1);
    }
  }
}

double PostInverseProblem::compFlowRate(std::vector<std::vector<double>> &vel)
{
  double flowRate;

  for (int k = 0; k < grid.nz + 1; k++)
  {
    for (int j = 0; j < grid.ny + 1; j++)
    {
      for (int i = 0; i < grid.nx + 1; i++)
      {
        int n = i + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
        if (crossSection == CrossSection::YZ)
        {
          if (i == crossPoint)
          {
            double area = grid.dy * grid.dz;
            flowRate += vel[n][flowRateVelDir] * area;
          }
        }
      }
    }
  }

  return flowRate;
}

double PostInverseProblem::compFlowRateError(const double flowRateRefVec, const double flowRateOptVec)
{
  double error;
  error = fabs(flowRateRefVec - flowRateOptVec) / flowRateRefVec;
  error *= 100;
  return error;
}