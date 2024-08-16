/**
 * @file GridCreation.h
 * @brief Header file for GridCreation class
 * @date August, 2024
 */

#ifndef GRIDCREATION_H
#define GRIDCREATION_H

#include <sys/stat.h>
#include <map>
#include <string>
#include <vector>
#include <memory>
#include "Config.h"
#include "Grid.h"
#include "Export.h"
#include "metis.h"

class GridCreation
{
public:
  GridCreation(Config &conf);

  void initialize(Config &conf);
  void divideWholeGrid();
  void collectLocalGrid();
  void outputDat();
  void outputVTU();

private:
  Grid grid;
  std::string outputDir;
  int CBExtraction;

  std::vector<int> cellId;
  std::vector<int> nodeId;
  std::vector<double> phi;

  std::map<int, std::vector<double>> vDirichlet;
  std::map<int, double> pDirichlet;

  std::vector<int> CBNodeMap;
  std::vector<int> CBCellMap;
  std::vector<std::vector<int>> CBNodeMapInCell;

  std::set<int> fluidUniqueNodes;

  void initializeCells(Config &conf);
  void initializeNodes(Config &conf);
  void setCellType(int nNodesInCell);
  void partitionMesh(int *eptr, int *eind, int nparts);
};

#endif // GRIDCREATION_H