/**
 * @file GridCreation.h
 * @author K.Ueda
 * @date August, 2024
 */

#include <sys/stat.h>
#include <map>
#include "Config.h"
#include "Grid.h"
#include "Export.h"

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

    std::string dir;
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/vtu";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  }
  
  Grid grid;

  std::string outputDir;
  int extractCB;

  std::vector<int> cellId;
  std::vector<int> nodeId;
  std::vector<double> phi;

  std::map<int, std::vector<double>> vDirichlet;
  std::map<int, double> pDirichlet;
  
  std::vector<int> mapCB;
  std::vector<int> mapCBCell;
  std::vector<std::vector<int>> mapCBInCell;
  
  void initialize(Config &conf);
  void divideWholeGrid();
  void collectLocalGrid();
  void outputDat();
  void outputVTU();
};