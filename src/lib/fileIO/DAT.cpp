/**
 * @file DAT.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "FileIO.h"

void DAT::exportCellDataDAT(const std::string &file, Cell &cell)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      ofs << cell(ic).node[p] << " ";
    }
    ofs << "\n";
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

void DAT::exportNodeDataDAT(const std::string &file, Node &node)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    for (int d = 0; d < 3; d++)
    {
      ofs << node.x[in][d] << " ";
    }
    ofs << "\n";
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}