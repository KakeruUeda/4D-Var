/**
 * @file GridFiltering.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/*****************************************
 * @brief Set velocity dirichlet boundary
 *        condition zero on solid nodes.
 */
void Config::setSolidDirichletValue()
{
  if(gridType != GridType::STRUCTURED)
  {
    return;
  }

  std::vector<double> vecTmp;
  vecTmp.resize(3, 0e0);

  for (int ic = 0; ic < nCellsGlobal; ic++)
  {
    if (phi[ic] < 1e-12)
    {
      for (int p = 0; p < nNodesInCell; p++)
      {
        std::vector<double> vecTmp(dim, 0e0);
        vDirichlet[cell[ic][p]] = vecTmp;
      }
    }
  }
}

void Config::getUniqueCells()
{
  for (int ic = 0; ic < phi.size(); ic++)
  {
    if (phi[ic] > 1e-12)
    {
      uniqueCells.insert(ic);
    }
  }
}

void Config::getUniqueNodes()
{
  for (int ic : uniqueCells)
  {
    for (int node : cell[ic])
    {
      uniqueNodes.insert(node);
    }
  }
}

void Config::getNewFilterdMap()
{
  int newIdx = 0;
  for (int oldIdx : uniqueCells)
  {
    cellMapping[oldIdx] = newIdx++;
  }

  newIdx = 0;
  for (int oldIdx : uniqueNodes)
  {
    nodeMapping[oldIdx] = newIdx++;
  }
}

void Config::filterCell()
{
  std::vector<std::vector<int>> filteredValues;
  filteredValues.reserve(uniqueCells.size());

  for (int ic : uniqueCells)
  {
    filteredValues.push_back(std::move(cell[ic]));
  }

  cell = std::move(filteredValues);
}

void Config::filterPhi()
{
  std::vector<double> filteredValues;
  filteredValues.reserve(uniqueCells.size());

  for (int ic : uniqueCells)
  {
    filteredValues.push_back(std::move(phi[ic]));
  }

  phi = std::move(filteredValues);
}

void Config::filterNode()
{
  std::vector<std::vector<double>> filterdValues;
  filterdValues.reserve(uniqueNodes.size());

  for (int in : uniqueNodes)
  {
    filterdValues.push_back(node[in]);
  }

  node = std::move(filterdValues);
}

void Config::filterVelocityDirichlet()
{
  std::map<int, std::vector<double>> fileterdValues;

  for (const auto &entry : vDirichlet)
  {
    if (uniqueNodes.find(entry.first) != uniqueNodes.end())
    {
      fileterdValues[entry.first] = entry.second;
    }
  }
  vDirichlet = std::move(fileterdValues);
}

void Config::filterPressureDirichlet()
{
  std::map<int, double> fileterdValues;

  for (const auto &entry : pDirichlet)
  {
    if (uniqueNodes.find(entry.first) != uniqueNodes.end())
    {
      fileterdValues[entry.first] = entry.second;
    }
  }
  pDirichlet = std::move(fileterdValues);
}

/***************************************************
 * @brief Extract fluid domain from structured grid
 *        to lower the computational cost.
 *        Might need to simplify.
 */
void Config::filterFluidGrid()
{
  if(!((gridType == GridType::STRUCTURED) && (extractFluid == ON)))
  {
    return;
  }

  getUniqueCells();
  getUniqueNodes();
  getNewFilterdMap();

  nCellsGlobal = uniqueCells.size();
  nNodesGlobal = uniqueNodes.size();

  filterCell();
  filterPhi();
  filterNode();

  filterVelocityDirichlet();
  filterPressureDirichlet();

  // New Cell Mapping
  for (auto &vec : cell)
  {
    for (auto &val : vec)
    {
      val = nodeMapping[val];
    }
  }

  // New vDirichlet Mapping
  std::map<int, std::vector<double>> vtmp;
  for (auto &entry : vDirichlet)
  {
    vtmp[nodeMapping[entry.first]] 
    = std::move(entry.second);
  }
  vDirichlet = std::move(vtmp);

  // New pDirichlet Mapping
  std::map<int, double> ptmp;
  for (auto &entry : pDirichlet)
  {
    ptmp[nodeMapping[entry.first]] 
    = std::move(entry.second);
  }
  pDirichlet = std::move(ptmp);

}
