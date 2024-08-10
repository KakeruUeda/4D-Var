/**
 * @file GridCreation.h
 * @author K.Ueda
 * @date August, 2024
 */

#include <sys/stat.h>
#include <map>
#include "Config.h"
#include "Grid.h"
#include "FileIO.h"

class GridCreation
{
public:
  GridCreation(Config &conf)
  {
    outputDir = conf.outputDir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + conf.outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  }
  Grid grid;

  std::string outputDir;

  std::vector<int> cellId;
  std::vector<int> nodeId;

  std::map<int, std::vector<double>> vDirichlet;
  std::map<int, double> pDirichlet;
  
  void initialize(Config &conf);
  void divideWholeGrid();
  void collectLocalGrid();
  void outputDat();
  void outputVTU();
};