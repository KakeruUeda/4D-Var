#include "Config.h"

/**
 * @brief Set velocity Dirichlet boundary
 *        condition to zero on solid nodes.
 */
void Config::setSolidDirichletValue()
{
  if(gridType != GridType::STRUCTURED) {
    return;
  }

  std::vector<double> vecTmp(dim, 0.0);

  for(int ic = 0; ic < nCellsGlobal; ic++) {
    if(phi[ic] < 1e-12) {
      for(int p = 0; p < nNodesInCell; p++) {
        vDirichlet[cell[ic][p]] = vecTmp;
      }
    }
  }
}

/**
 * @brief Identify unique cells.
 */
void Config::getUniqueCells()
{
  for(int ic = 0; ic < phi.size(); ic++) {
    if(phi[ic] > 0) {
      fluidUniqueCells.insert(ic);
    }
  }
}

/**
 * @brief Identify unique nodes.
 */
void Config::getUniqueNodes()
{
  for(int ic : fluidUniqueCells) {
    for(int node : cell[ic]) {
      fluidUniqueNodes.insert(node);
    }
  }
}

/**
 * @brief Identify unique control boundary cells.
 */
void Config::getUniqueCBCells()
{
  for(int cb : CBCellMap) {
    if(phi[cb] > 1e-12) {
      fluidUniqueCBCells.insert(cb);
    }
  }
}

/**
 * @brief Identify unique control boundary cells index.
 */
void Config::getUniqueCBCellsIdx()
{
  for(int icb = 0; icb < CBCellMap.size(); icb++) {
    int cb = CBCellMap[icb];
    if(phi[cb] > 1e-12) {
      fluidUniqueCBCellsIdx.insert(icb);
    }
  }
}

/**
 * @brief Identify unique control boundary nodes.
 */
void Config::getUniqueCBNodes()
{
  for(int icb = 0; icb < CBCellMap.size(); icb++) {
    int cb = CBCellMap[icb];
    if(phi[cb] > 1e-12) {
      for(int node : CBNodeMapInCell[icb]) {
        fluidUniqueCBNodes.insert(node);
      }
    }
  }

  std::set<int> tmp;
  tmp = fluidUniqueCBNodes;

  for(int icb = 0; icb < CBCellMap.size(); icb++) {
    int cb = CBCellMap[icb];
    if(phi[cb] < 1e-12) {
      for(int node : CBNodeMapInCell[icb]) {
        int erased = tmp.erase(node);
        if(erased > 0) {
          fluidUniqueCBEdgeNodes.insert(node);
        }
      }
    }
  }
}

/**
 * @brief Create new filtered map for cells and nodes.
 */
void Config::getNewFilterdMap()
{
  int newIdx = 0;
  for(int oldIdx : fluidUniqueCells) {
    cellMapping[oldIdx] = newIdx++;
  }

  newIdx = 0;
  for(int oldIdx : fluidUniqueNodes) {
    nodeMapping[oldIdx] = newIdx++;
  }
}

/** 
 * @brief Filter cells based on unique cells.
 */
void Config::filterCell()
{
  std::vector<std::vector<int>> filteredValues;

  for(int ic : fluidUniqueCells) {
    filteredValues.push_back(cell[ic]);
  }

  cell = std::move(filteredValues);
}

/**
 * @brief Filter phi values based on unique cells.
 */
void Config::filterPhi()
{
  std::vector<double> filteredValues;

  for(int ic : fluidUniqueCells) {
    filteredValues.push_back(phi[ic]);
  }

  phi = std::move(filteredValues);
}

/**
 * @brief Filter nodes based on unique nodes.
 */
void Config::filterNode()
{
  std::vector<std::vector<double>> filteredValues;

  for(int in : fluidUniqueNodes) {
    filteredValues.push_back(node[in]);
  }

  node = std::move(filteredValues);
}

/**
 * @brief Filter control boundary nodes.
 */
void Config::filterMapCBNode()
{
  std::vector<int> filteredValues;

  for(int in : fluidUniqueCBNodes) {
    filteredValues.push_back(in);
  }

  CBNodeMap = std::move(filteredValues);
}

/**
 * @brief Filter nodes based on unique nodes.
 */
void Config::filterMapCBEdgeNode()
{
  std::vector<int> filteredValues;

  for(int in : fluidUniqueCBEdgeNodes) {
    filteredValues.push_back(in);
  }

  CBEdgeNodeMap = std::move(filteredValues);
}

/**
 * @brief Filter control boundary cells.
 */
void Config::filterMapCBCell()
{
  std::vector<int> filteredValues;

  for(int cb : fluidUniqueCBCells) {
    filteredValues.push_back(cb);
  }

  CBCellMap = std::move(filteredValues);
}

/**
 * @brief Filter control boundary nodes in cells.
 */
void Config::filterMapCBInCell()
{
  std::vector<std::vector<int>> filteredValues;
  for(int cb : fluidUniqueCBCellsIdx) {
    filteredValues.push_back(CBNodeMapInCell[cb]);
  }

  CBNodeMapInCell = std::move(filteredValues);
}

/**
 * @brief Filter velocity Dirichlet boundary conditions.
 */
void Config::filterVelocityDirichlet()
{
  std::map<int, std::vector<double>> filteredValues;

  for(const auto &entry : vDirichlet) {
    if(fluidUniqueNodes.find(entry.first) != fluidUniqueNodes.end()) {
      filteredValues[entry.first] = entry.second;
    }
  }
  vDirichlet = std::move(filteredValues);
}

/**
 * @brief Filter pressure Dirichlet boundary conditions.
 */
void Config::filterPressureDirichlet()
{
  std::map<int, double> filteredValues;

  for(const auto &entry : pDirichlet) {
    if(fluidUniqueNodes.find(entry.first) != fluidUniqueNodes.end()) {
      filteredValues[entry.first] = entry.second;
    }
  }
  pDirichlet = std::move(filteredValues);
}

/**
 * @brief Apply new cell and node mappings.
 */
void Config::applyMapping()
{
  for(auto &vec : cell) {
    for(auto &value : vec) {
      value = nodeMapping[value];
    }
  }

  for(auto &value : CBCellMap) {
    value = cellMapping[value];
  }

  for(auto &vec : CBNodeMapInCell) {
    for(auto &value : vec) {
      value = nodeMapping[value];
    }
  }

  for(auto &value : CBNodeMap) {
    value = nodeMapping[value];
  }

  for(auto &value : CBEdgeNodeMap) {
    value = nodeMapping[value];
  }

  std::map<int, std::vector<double>> vtmp;
  for(auto &entry : vDirichlet) {
    vtmp[nodeMapping[entry.first]] = std::move(entry.second);
  }
  vDirichlet = std::move(vtmp);

  std::map<int, double> ptmp;
  for(auto &entry : pDirichlet) {
    ptmp[nodeMapping[entry.first]] = std::move(entry.second);
  }
  pDirichlet = std::move(ptmp);
}

/**
 * @brief Filter fluid domain from structured grid.
 */
void Config::filterFluidGrid()
{
  if(!((gridType == GridType::STRUCTURED) && (fluidExtraction == ON))) {
    return;
  }
  // Get unique cells and nodes
  getUniqueCells();
  getUniqueNodes();

  if(CBExtraction == ON) {
    getUniqueCBCells();
    getUniqueCBNodes();
    getUniqueCBCellsIdx();
  }

  getNewFilterdMap();

  nCellsGlobal = fluidUniqueCells.size();
  nNodesGlobal = fluidUniqueNodes.size();

  // Filtering
  filterCell();
  filterPhi();
  filterNode();

  if(CBExtraction == ON) {
    filterMapCBNode();
    filterMapCBEdgeNode();
    filterMapCBCell();
    filterMapCBInCell();
  }

  filterVelocityDirichlet();
  filterPressureDirichlet();

  // Apply new mappings (node numbers start from 0)
  applyMapping();
}
