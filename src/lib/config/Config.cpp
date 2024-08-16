/**
 * @file Config.cpp
 * @author K.Ueda
 * @date May, 2024
 */

#include "Config.h"

/**************************
 * @brief Construct config.
 */
Config::Config(std::string inputFile, std::string appName)
{
  setApplication(appName);
  tryOpenConfigFile(inputFile);
  tryReadConfigFile();
}

/**************************
 * @brief Set application.
 */
void Config::setApplication(std::string appName)
{
  if(appName == "SNS") {
    app = Application::SNS;
  } else if(appName == "USNS") {
    app = Application::USNS;
  } else if(appName == "VOXELDATACREATION") {
    app = Application::VOXELDATACREATION;
  } else if(appName == "FDVAR") {
    app = Application::FDVAR;
  } else if(appName == "FLOWRATE") {
    app = Application::FLOWRATE;
  } else if(appName == "MAE") {
    app = Application::MAE;
  } else if(appName == "GRIDCREATION") {
    app = Application::GRIDCREATION;
  } else if(mpi.myId == 0) {
    std::cout << "Unknown appName" << std::endl;
  }
}

/**************************
 * @brief Open config file.
 */
void Config::tryOpenConfigFile(std::string inputFile)
{
  try {
    int error;
    if((error = tp.read(inputFile)) != TP_NO_ERROR) {
      throw std::runtime_error("Open error");
    }
  } catch(const std::runtime_error &e) {
    if(mpi.myId == 0) {
      std::cout << e.what() << std::endl;
    }
    isReadingError = true;
  }
}

/**************************
 * @brief Read config file.
 */
void Config::tryReadConfigFile()
{
  try {
    readConfigFile();
  } catch(const std::runtime_error &e) {
    if(mpi.myId == 0) {
      std::cout << e.what() << std::endl;
    }
    isReadingError = true;
  }
}

/**************************
 * @brief Read config file.
 */
void Config::readConfigFile()
{
  switch(app) {
  case Application::GRIDCREATION: {
    TextReaderGridCreation reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readStructuredBoundaryInfo(*this);
    break;
  }
  case Application::USNS: {
    TextReaderUSNS reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readBoundaryInfo(*this);
    reader.readPhysicalInfo(*this);
    reader.readDarcyInfo(*this);
    reader.readTimeInfo(*this);
    break;
  }
  case Application::VOXELDATACREATION: {
    TextReaderVoxelDataCreation reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readSnapInfo(*this);
    reader.readOriginalInfo(*this);
    break;
  }
  case Application::FDVAR: {
    TextReader4DVar reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readBoundaryInfo(*this);
    reader.readPhysicalInfo(*this);
    reader.readDarcyInfo(*this);
    reader.readTimeInfo(*this);
    reader.readInverseInfo(*this);
    reader.readDataInfo(*this);
    break;
  }

  case Application::FLOWRATE:
    readBasicParameter();
    readPostInverseBasicParameter();
    readPostInverseVelocityParameter();
    readPostInverseFlowRateParameter();
    break;

  default:
    throw std::runtime_error("Unknown Application");
    break;
  }
}

/*****************************************
 * @brief Set velocity dirichlet boundary
 *        condition zero on solid nodes.
 */
void Config::setSolidBoundary()
{
  std::vector<double> vecTmp;
  vecTmp.resize(3, 0e0);

  for(int ic = 0; ic < nCellsGlobal; ic++) {
    if(phi[ic] < 1e-12) {
      for(int p = 0; p < nNodesInCell; p++) {
        std::vector<double> vecTmp(dim, 0e0);
        vDirichlet[cell[ic][p]] = vecTmp;
      }
    }
  }

  bool flag;
  for(int ic = 0; ic < nCellsGlobal; ic++) {
    flag = false;
    if(phi[ic] < 1e-12) {
      for(int ic2 = 0; ic2 < controlCellMap.size(); ic2++) {
        if(ic == controlCellMap[ic2]) {
          flag = true;
        }
      }
      if(flag == false) {
        for(int p = 0; p < nNodesInCell; p++) {
          std::vector<double> vecTmp(dim, 0e0);
          vDirichletWall[cell[ic][p]] = vecTmp;
        }
      }
    }
  }
}

/***************************************************
 * @brief Extract fluid domain from structured grid
 *        to lower computational cost.
 *        Might need to simplify.
 */
void Config::setFluidDomain()
{
  int nCellsGlobalTmp = nCellsGlobal;
  std::vector<std::vector<int>> cellTmp = cell;

  nCellsGlobal = 0;
  cell.clear();
  for(int ic = 0; ic < nCellsGlobalTmp; ic++) {
    if(phi[ic] < 1e-12) {
      continue;
    }
    sortCell.push_back(ic);
    cell.push_back(cellTmp[ic]);
    nCellsGlobal++;
  }

  std::vector<std::vector<int>> controlNodeInCellTmp = controlNodeInCell;
  controlNodeInCell.clear();
  for(int ic = 0; ic < controlNodeInCellTmp.size(); ic++) {
    if(phi[controlCellMap[ic]] < 1e-12) {
      continue;
    }
    controlNodeInCell.push_back(controlNodeInCellTmp[ic]);
  }

  std::vector<int> controlCellMapTmp = controlCellMap;
  controlCellMap.clear();
  for(int ic = 0; ic < controlCellMapTmp.size(); ic++) {
    if(phi[controlCellMapTmp[ic]] < 1e-12) {
      continue;
    }
    controlCellMap.push_back(controlCellMapTmp[ic]);
  }

  int count2 = 0;
  std::vector<int> convertCellOldToNew(nCellsGlobalTmp, 0);
  for(int in = 0; in < sortCell.size(); in++) {
    convertCellOldToNew[sortCell[in]] = count2;
    count2++;
  }

  for(int ic = 0; ic < controlCellMap.size(); ic++) {
    controlCellMap[ic] = convertCellOldToNew[controlCellMap[ic]];
  }

  std::vector<double> phiTmp = phi;
  phi.clear();
  int count = 0;
  for(int ic = 0; ic < nCellsGlobalTmp; ic++) {
    if(phiTmp[ic] < 1e-12) {
      continue;
    }
    phi.push_back(phiTmp[ic]);
  }

  sortNode.resize(nCellsGlobal * nNodesInCell);

  count = 0;
  for(int ic = 0; ic < nCellsGlobal; ic++) {
    for(int p = 0; p < nNodesInCell; p++) {
      sortNode[count++] = cell[ic][p];
    }
  }

  sort(sortNode.begin(), sortNode.end());
  sortNode.erase(unique(sortNode.begin(), sortNode.end()), sortNode.end());

  int nNodesGlobalTmp = nNodesGlobal;
  nNodesGlobal = 0;
  std::vector<int> sortNodeNew(sortNode.size(), 0);
  std::vector<int> convertNodeOldToNew(nNodesGlobalTmp, 0);
  for(int in = 0; in < sortNode.size(); in++) {
    sortNodeNew[in] = nNodesGlobal++;
    convertNodeOldToNew[sortNode[in]] = sortNodeNew[in];
  }

  for(int ic = 0; ic < nCellsGlobal; ic++) {
    for(int p = 0; p < nNodesInCell; p++) {
      cell[ic][p] = convertNodeOldToNew[cell[ic][p]];
    }
  }

  for(int ic = 0; ic < controlNodeInCell.size(); ic++) {
    for(int p = 0; p < nControlNodesInCell; p++) {
      controlNodeInCell[ic][p] = convertNodeOldToNew[controlNodeInCell[ic][p]];
    }
  }

  std::vector<std::vector<double>> nodeTmp = node;
  node.erase(node.begin(), node.end());

  for(int in = 0; in < sortNode.size(); in++) {
    node.push_back(nodeTmp[sortNode[in]]);
  }

  std::map<int, std::vector<double>> vDirichletTmp = vDirichlet;
  std::map<int, std::vector<double>> vDirichletWallTmp = vDirichletWall;
  std::map<int, double> pDirichletTmp = pDirichlet;

  vDirichlet.clear();
  vDirichletWall.clear();
  pDirichlet.clear();

  for(int in = 0; in < sortNode.size(); in++) {
    if(vDirichletTmp.size() == 0) {
      break;
    }
    int key = sortNode[in];
    auto it = vDirichletTmp.find(key);
    if(it != vDirichletTmp.end()) {
      int keyNew = convertNodeOldToNew[key];
      vDirichlet[keyNew] = vDirichletTmp[key];
    }
  }

  for(int in = 0; in < sortNode.size(); in++) {
    if(vDirichletWallTmp.size() == 0) {
      break;
    }
    int key = sortNode[in];
    auto it = vDirichletWallTmp.find(key);
    if(it != vDirichletWallTmp.end()) {
      int keyNew = convertNodeOldToNew[key];
      vDirichletWall[keyNew] = vDirichletWallTmp[key];
    }
  }

  for(int in = 0; in < sortNode.size(); in++) {
    if(pDirichletTmp.size() == 0) {
      break;
    }
    int key = sortNode[in];
    auto it = pDirichletTmp.find(key);
    if(it != pDirichletTmp.end()) {
      int keyNew = convertNodeOldToNew[key];
      pDirichlet[keyNew] = pDirichletTmp[key];
    }
  }

  std::vector<int> controlBoundaryMapTmp = controlBoundaryMap;
  std::vector<int> controlBoundaryMapTmp2;
  controlBoundaryMap.clear();

  for(int in = 0; in < sortNode.size(); in++) {
    if(controlBoundaryMapTmp.size() == 0) {
      break;
    }
    int key = sortNode[in];
    if(std::find(controlBoundaryMapTmp.begin(), controlBoundaryMapTmp.end(), key) != controlBoundaryMapTmp.end()) {
      controlBoundaryMap.push_back(convertNodeOldToNew[key]);
      controlBoundaryMapTmp2.push_back(key);
    }
  }

  /*
		 for(int ib=0; ib<controlBoundaryMapTmp.size(); ib++){
			bool flag = true;
			bool flag2 = false;
			for(int ic=0; ic<nCellsGlobalTmp; ic++){
					if(phiTmp[ic] < 1e-12){
							for(int p=0; p<nNodesInCell; p++){
									if(controlBoundaryMapTmp[ib] == cellTmp[ic][p]){
											flag = false;
											flag2 = true;
									}
							}
					}
					if(flag2 == true) break;
			}
			if(flag == true){
					int n = controlBoundaryMapTmp[ib];
					controlBoundaryMap.push_back(convertNodeOldToNew[n]);
			}
		 }
	 */
  int count22 = 0;

  isBoundaryEdge.resize(nNodesGlobal, false);
  for(int ic = 0; ic < cellTmp.size(); ic++) {
    if(phiTmp[ic] < 1e-12) {
      for(int p = 0; p < nNodesInCell; p++) {
        int key = cellTmp[ic][p];
        if(std::find(controlBoundaryMapTmp2.begin(), controlBoundaryMapTmp2.end(), key) !=
           controlBoundaryMapTmp2.end()) {
          isBoundaryEdge[convertNodeOldToNew[key]] = true;
        }
      }
    }
  }
}
