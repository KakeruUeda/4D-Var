/**
 * @file ReadTP.cpp
 * @author K.Ueda
 * @date May, 2024
 */

#include "Config.h"

/*****************************
 * @brief Read text parameter.
 */
void Config::readBasicParameter()
{
  std::string str, base_label, label;

  base_label = "/Base";
  label = base_label + "/dim";

  if(!tp.getInspectedValue(label, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/outputDir";
  if(!tp.getInspectedValue(label, outputDir)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/numOfOMP";
  if(!tp.getInspectedValue(label, nOMP)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readPysicalParameter()
{
  std::string str, base_label, label;

  base_label = "/PysicalParameter";
  label = base_label + "/rho";
  if(!tp.getInspectedValue(label, rho)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/mu";
  if(!tp.getInspectedValue(label, mu)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readNRParameter()
{
  std::string str, base_label, label;

  base_label = "/NRParameter";
  label = base_label + "/NRtolerance";
  if(!tp.getInspectedValue(label, NRtolerance)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readTimeParameter()
{
  std::string str, base_label, label;

  base_label = "/TimeParameter";
  label = base_label + "/dt";
  if(!tp.getInspectedValue(label, dt)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/timeMax";
  if(!tp.getInspectedValue(label, timeMax)) {
    throw std::runtime_error(label + " is not set");
  }

  std::string ON_OFF;

  label = base_label + "/pulsatileFlow";
  if(!tp.getInspectedValue(label, ON_OFF)) {
    throw std::runtime_error(label + " is not set");
  }

  if(ON_OFF == "ON") {
    pulsatileFlow = ON;
  } else if(ON_OFF == "OFF") {
    pulsatileFlow = OFF;
  } else {
    throw std::runtime_error("ON or OFF is not set");
  }

  label = base_label + "/pulseBeginItr";
  if(!tp.getInspectedValue(label, pulseBeginItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/T";
  if(!tp.getInspectedValue(label, T)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readDarcyParameter()
{
  std::string str, base_label, label;
  base_label = "/DarcyParameter";

  label = base_label + "/alpha";
  if(!tp.getInspectedValue(label, alpha)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/resistance";
  if(!tp.getInspectedValue(label, resistance)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readGridType()
{
  std::string str, base_label, label;
  std::string nodeFile;

  base_label = "/Grid";

  label = base_label + "/type";

  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str != "Structured" && str != "Unstructured") {
    throw std::runtime_error("Unknown GridType");
  }

  if(str == "Structured") {
    gridType = GridType::STRUCTURED;
  } else if(str == "Unstructured") {
    gridType = GridType::UNSTRUCTURED;
  }

  label = base_label + "/nNodesInCell";
  if(!tp.getInspectedValue(label, nNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readGridParameter()
{
  std::string str, base_label, sub_label, label;

  base_label = "/Grid";
  sub_label = base_label + "/GridType";

  label = sub_label + "/type";
  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "Structured") {
    gridType = GridType::STRUCTURED;
  } else if(str == "Unstructured") {
    gridType = GridType::UNSTRUCTURED;
  } else {
		throw std::runtime_error("Unknown GridType");
	}

  label = sub_label + "/nNodesInCell";
  if(!tp.getInspectedValue(label, nNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }

  switch(app) {
  case Application::GRIDCREATION:
    readStrGridParameter();
    break;
  case Application::USNS:
    readBaseGridParameter();
    readSubGridParameter();
    if(gridType == GridType::STRUCTURED) {
      readStrGridParameter();
    }
  default:
    throw std::runtime_error("Application is not set");
  }
}

void Config::readBaseGridParameter()
{
  std::string str, base_label, sub_label, label;

  base_label = "/Grid";
  sub_label = base_label + "/BaseGrid";
  label = sub_label + "/node";

  std::string nodeFile;
  if(!tp.getInspectedValue(label, nodeFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorDataDAT<double>(nodeFile, node);
  nNodesGlobal = node.size();

  std::string cellFile;
  label = sub_label + "/cell";
  if(!tp.getInspectedValue(label, cellFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorDataDAT<int>(cellFile, cell);
  nCellsGlobal = cell.size();
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readSubGridParameter()
{
  std::string str, base_label, sub_label, label;
  std::string nodeIdFile;
  std::string cellIdFile;

  base_label = "/Grid";

  sub_label = base_label + "/SubGrid";

  label = sub_label + "/nodeId";
  if(!tp.getInspectedValue(label, nodeIdFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<int>(nodeIdFile, nodeId);

  label = sub_label + "/cellId";
  if(!tp.getInspectedValue(label, cellIdFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<int>(cellIdFile, cellId);
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readStrGridParameter()
{
  std::string str, base_label, sub_label, label;
  int tmpInt[dim];
  double tmpDouble[dim];

  base_label = "/Grid";

  sub_label = base_label + "/StructuredGrid";
  label = sub_label + "/nx";

  if(!tp.getInspectedVector(label, tmpInt, dim)) {
    throw std::runtime_error(label + " is not set.");
  }

  nx = tmpInt[0];
  ny = tmpInt[1];
  nz = tmpInt[2];

  label = sub_label + "/lx";
  if(!tp.getInspectedVector(label, tmpDouble, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  lx = tmpDouble[0];
  ly = tmpDouble[1];
  lz = tmpDouble[2];

  dx = lx / (double)nx;
  dy = ly / (double)ny;
  dz = lz / (double)nz;

  label = sub_label + "/image";

  std::string imageFile;
  if(!tp.getInspectedValue(label, imageFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<double>(imageFile, phi);

  label = sub_label + "/fluidExtraction";
  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(app == Application::GRIDCREATION) {
    str = "OFF";
  }

  if(str == "ON") {
    fluidExtraction = ON;
  } else if(str == "OFF") {
    fluidExtraction = OFF;
  } else {
    throw std::runtime_error("ON or OFF is not set");
  }

  if(fluidExtraction) {
    label = sub_label + "/fluidUniqueNodes";
    std::string fluidGridFile;
    if(!tp.getInspectedValue(label, fluidGridFile)) {
      throw std::runtime_error(label + " is not set");
    }
    IMPORT::importScalarDataDAT<int>(fluidGridFile, vecFluidUniqueNodes);
  }else{
		vecFluidUniqueNodes.resize(0);
	}

  if(app == Application::GRIDCREATION) {
    setStrGrid(); // define str grid
  } 
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readStrBoundaryParameter()
{
  std::string str, base_label;
  std::string face;
  std::string labelFace, labelType, labelValue;
  std::string bdTypeTmp;
  int tmp = 0;

  base_label = "/StructuredBoundary";

  face = "bottom";
  labelFace = base_label + "/bottom";
  labelType = base_label + "/bottom/type";
  labelValue = base_label + "/bottom/value";
  readStrBoundaryValue(face, labelFace, labelType, labelValue);

  face = "top";
  labelFace = base_label + "/top";
  labelType = base_label + "/top/type";
  labelValue = base_label + "/top/value";
  readStrBoundaryValue(face, labelFace, labelType, labelValue);

  face = "left";
  labelFace = base_label + "/left";
  labelType = base_label + "/left/type";
  labelValue = base_label + "/left/value";
  readStrBoundaryValue(face, labelFace, labelType, labelValue);

  face = "right";
  labelFace = base_label + "/right";
  labelType = base_label + "/right/type";
  labelValue = base_label + "/right/value";
  readStrBoundaryValue(face, labelFace, labelType, labelValue);

  face = "front";
  labelFace = base_label + "/front";
  labelType = base_label + "/front/type";
  labelValue = base_label + "/front/value";
  readStrBoundaryValue(face, labelFace, labelType, labelValue);

  face = "back";
  labelFace = base_label + "/back";
  labelType = base_label + "/back/type";
  labelValue = base_label + "/back/value";
  readStrBoundaryValue(face, labelFace, labelType, labelValue);
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readStrBoundaryValue(std::string face, std::string labelFace, std::string labelType,
                                  std::string labelValue)
{
  std::string bdTypeTmp;
  if(!tp.getInspectedValue(labelType, bdTypeTmp)) {
    throw std::runtime_error(labelType + " is not set");
  }

  bdType.push_back(bdTypeTmp);

  if(bdTypeTmp == "v") {
    double value[dim];
    if(!tp.getInspectedVector(labelValue, value, dim)) {
      throw std::runtime_error(labelValue + " is not set");
    }

    setBoundaryVelocityValue(face, value);
  } else if(bdTypeTmp == "p") {
    double value;
    if(!tp.getInspectedValue(labelValue, value)) {
      throw std::runtime_error(labelValue + " is not set");
    }

    setBoundaryPressureValue(face, value);
  } else if(bdTypeTmp == "poiseuille") {
    std::string label = labelFace + "/center";

    if(!tp.getInspectedVector(label, center, 3)) {
      throw std::runtime_error(label + " is not set");
    }

    label = labelFace + "/R";

    if(!tp.getInspectedValue(label, R)) {
      throw std::runtime_error(label + " is not set");
    }

    label = labelFace + "/Q";

    if(!tp.getInspectedValue(label, Q)) {
      throw std::runtime_error(label + " is not set");
    }

    setBoundaryPoiseuilleValue(face);
  } else if(bdTypeTmp == "free") {
  } else {
    throw std::runtime_error("label " + bdTypeTmp + " undefined");
  }

  return;
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readBoundaryParameter()
{
  std::string str, base_label, label;
  std::string nodeFile;

  base_label = "/Boundary";

  std::string velFile;
  label = base_label + "/velocityDirichlet";

  if(!tp.getInspectedValue(label, velFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::ifstream ifsVel(velFile);
  if(!ifsVel.is_open()) {
    throw std::runtime_error("Failed to open file: " + velFile);
  }
  while(getline(ifsVel, str)) {
    int index;
    std::istringstream iss(str);
    std::vector<double> vecTmp;

    for(int d = 0; d < dim + 1; d++) {
      getline(iss, str, ' ');
      if(d == 0) {
        index = stoi(str);
      } else {
        vecTmp.push_back(stod(str));
      }
    }
    vDirichlet[index] = vecTmp;
  }
  ifsVel.close();

  std::string preFile;
  label = base_label + "/pressureDirichlet";

  if(!tp.getInspectedValue(label, preFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::ifstream ifsPre(preFile);
  if(!ifsPre.is_open()) {
    throw std::runtime_error("Failed to open file: " + preFile);
  }
  while(getline(ifsPre, str)) {
    int index;
    std::istringstream iss(str);
    std::vector<double> preTmp;

    for(int d = 0; d < 1 + 1; d++) {
      getline(iss, str, ' ');
      if(d == 0) {
        index = stoi(str);
      } else {
        pDirichlet[index] = stoi(str);
      }
    }
  }
  ifsPre.close();
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readControlBoundaryParameter()
{
  std::string str, base_label, label;
  base_label = "/StructuredBoundary";

  label = base_label + "/controlBoundaryExtraction";

  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "ON") {
    CBExtraction = ON;
  } else if(str == "OFF") {
    CBExtraction = OFF;
    return;
  } else {
    throw std::runtime_error("ON or OFF is not set");
  }

  label = base_label + "/inletControlBoundary";

  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "left") {
    inletCB = ControlBoundary::left;
  }
  if(str == "right") {
    inletCB = ControlBoundary::right;
  }
  if(str == "top") {
    inletCB = ControlBoundary::top;
  }
  if(str == "bottom") {
    inletCB = ControlBoundary::bottom;
  }
  if(str == "front") {
    inletCB = ControlBoundary::front;
  }
  if(str == "back") {
    inletCB = ControlBoundary::back;
  }

  setControlBoundary();
}


/*****************************
 * @brief Read text parameter.
 */
void Config::readInverseParameter()
{
  std::string str, base_label, label;
  base_label = "/Inverse";

  label = base_label + "/controlBoundary";

  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  planeDir.resize(2, 0);

  if(str == "left") {
    controlBoundary = ControlBoundary::left;
    planeDir[0] = 1, planeDir[1] = 2;
  } else if(str == "right") {
    controlBoundary = ControlBoundary::right;
    planeDir[0] = 1, planeDir[1] = 2;
  } else if(str == "top") {
    controlBoundary = ControlBoundary::top;
    planeDir[0] = 0, planeDir[1] = 2;
  } else if(str == "bottom") {
    controlBoundary = ControlBoundary::bottom;
    planeDir[0] = 0, planeDir[1] = 2;
  } else if(str == "front") {
    controlBoundary = ControlBoundary::front;
    planeDir[0] = 0, planeDir[1] = 1;
  } else if(str == "back") {
    controlBoundary = ControlBoundary::back;
    planeDir[0] = 0, planeDir[1] = 1;
  }

  label = base_label + "/voxelVelocity";

  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "average") {
    vvox = VoxelVelocity::AVERAGE;
  } else if(str == "interpolation") {
    vvox = VoxelVelocity::INTERPOLATION;
  } else {
    throw std::runtime_error("undefined voxelVelocity");
  }

  label = base_label + "/aCF";
  if(!tp.getInspectedValue(label, aCF)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/bCF";
  if(!tp.getInspectedValue(label, bCF)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/gCF";
  if(!tp.getInspectedValue(label, gCF)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/alphaX0";
  if(!tp.getInspectedValue(label, alphaX0)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/alphaX";
  if(!tp.getInspectedValue(label, alphaX)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/outputItr";
  if(!tp.getInspectedValue(label, outputItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/loopMax";
  if(!tp.getInspectedValue(label, loopMax)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readDataParameter()
{
  std::string str, base_label, label;
  base_label = "/Data";

  int tmpInt[dim];
  double tmpDouble[dim];

  label = base_label + "/nControlNodesInCell";
  if(!tp.getInspectedValue(label, nControlNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nSnapShot";
  if(!tp.getInspectedValue(label, nSnapShot)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/snapInterval";
  if(!tp.getInspectedValue(label, snapInterval)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/snapTimeBeginItr";
  if(!tp.getInspectedValue(label, snapTimeBeginItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nNodesInDataCell";
  if(!tp.getInspectedValue(label, nNodesInCellData)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nxData";
  if(!tp.getInspectedVector(label, tmpInt, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  nxData = tmpInt[0];
  nyData = tmpInt[1];
  nzData = tmpInt[2];

  label = base_label + "/lxData";
  if(!tp.getInspectedVector(label, tmpDouble, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  lxData = tmpDouble[0];
  lyData = tmpDouble[1];
  lzData = tmpDouble[2];

  dxData = lxData / (double)nxData;
  dyData = lyData / (double)nyData;
  dzData = lzData / (double)nzData;

  nDataCellsGlobal = nxData * nyData * nzData;

  std::string controlBoundaryFile;
  label = base_label + "/controlBoundary";

  if(!tp.getInspectedValue(label, controlBoundaryFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::ifstream ifsControlBoundary(controlBoundaryFile);
  if(!ifsControlBoundary.is_open()) {
    throw std::runtime_error("Failed to open file: " + controlBoundaryFile);
  }
  while(getline(ifsControlBoundary, str)) {
    std::istringstream iss(str);
    getline(iss, str, ' ');
    controlBoundaryMap.push_back(stoi(str));
  }
  ifsControlBoundary.close();

  std::string controlCellMapFile;
  label = base_label + "/controlCellMap";

  if(!tp.getInspectedValue(label, controlCellMapFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::ifstream ifscontrolCellMap(controlCellMapFile);
  if(!ifscontrolCellMap.is_open()) {
    throw std::runtime_error("Failed to open file: " + controlCellMapFile);
  }
  while(getline(ifscontrolCellMap, str)) {
    std::istringstream iss(str);
    getline(iss, str, ' ');
    controlCellMap.push_back(stoi(str));
  }
  ifscontrolCellMap.close();

  std::string controlNodeInCellFile;
  label = base_label + "/controlNodeInCell";

  if(!tp.getInspectedValue(label, controlNodeInCellFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::ifstream ifsControlNodeInCell(controlNodeInCellFile);
  if(!ifsControlNodeInCell.is_open()) {
    throw std::runtime_error("Failed to open file: " + controlNodeInCellFile);
  }
  while(getline(ifsControlNodeInCell, str)) {
    std::istringstream iss(str);
    std::vector<int> ctrTmp;
    for(int p = 0; p < nControlNodesInCell; p++) {
      getline(iss, str, ' ');
      ctrTmp.push_back(stod(str));
    }
    controlNodeInCell.push_back(ctrTmp);
  }
  ifsControlNodeInCell.close();

  label = base_label + "/inputDir";
  if(!tp.getInspectedValue(label, inputDir)) {
    throw std::runtime_error(label + " is not set");
  }

  velocityData.resize(nSnapShot);

  // Read data from bin
  for(int step = 0; step < nSnapShot; step++) {
    std::string velFile = inputDir + "/data_" + std::to_string(step) + ".bin";
    IMPORT::importVectorDataBIN<double>(velFile, velocityData[step]);
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readPostInverseBasicParameter()
{
  std::string str, base_label, label;
  base_label = "/PostInverseBasic";

  label = base_label + "/gridType";
  std::string gridTypeString;

  if(!tp.getInspectedValue(label, gridTypeString)) {
    throw std::runtime_error(label + " is not set");
  }

  if(gridTypeString != "Structured" && gridTypeString != "Unstructured") {
    throw std::runtime_error("Unknown GridType");
  }

  if(gridTypeString == "Structured") {
    gridType = GridType::STRUCTURED;
  } else if(gridTypeString == "Unstructured") {
    gridType = GridType::UNSTRUCTURED;
  }

  if(gridType == GridType::STRUCTURED) {
    int tmpInt[dim];
    double tmpDouble[dim];

    label = base_label + "/nx";
    if(!tp.getInspectedVector(label, tmpInt, dim)) {
      throw std::runtime_error(label + " is not set.");
    }

    nx = tmpInt[0];
    ny = tmpInt[1];
    nz = tmpInt[2];

    label = base_label + "/lx";
    if(!tp.getInspectedVector(label, tmpDouble, dim)) {
      throw std::runtime_error(label + " is not set");
    }

    lx = tmpDouble[0];
    ly = tmpDouble[1];
    lz = tmpDouble[2];

    dx = lx / (double)nx;
    dy = ly / (double)ny;
    dz = lz / (double)nz;
  }

  label = base_label + "/nNodesInCell";
  if(!tp.getInspectedValue(label, nNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nRef";
  if(!tp.getInspectedValue(label, nRef)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readPostInverseVelocityParameter()
{
  std::string str, base_label, label;
  base_label = "/PostInverseVelocity";

  int num = 0;
  std::string velRefFile;
  velRef.resize(nRef);

  while(1) {
    if(num >= nRef) {
      break;
    }
    label = base_label + "/velRef" + std::to_string(num);

    if(!tp.getInspectedValue(label, velRefFile)) {
      throw std::runtime_error(label + " is not set");
    }

    std::ifstream ifsVelRef(velRefFile);
    if(!ifsVelRef.is_open()) {
      throw std::runtime_error("Failed to open file: " + velRefFile);
    }
    std::vector<std::vector<double>> velRefTmp;

    while(getline(ifsVelRef, str)) {
      std::istringstream iss(str);
      std::vector<double> tmp;
      for(int d = 0; d < dim; d++) {
        getline(iss, str, ' ');
        tmp.push_back(stod(str));
      }
      velRefTmp.push_back(tmp);
    }
    ifsVelRef.close();
    velRef[num] = velRefTmp;
    num++;
  }

  num = 0;
  std::string velOptFile;
  velOpt.resize(nRef);

  while(1) {
    if(num >= nRef) {
      break;
    }
    label = base_label + "/velOpt" + std::to_string(num);

    if(!tp.getInspectedValue(label, velOptFile)) {
      throw std::runtime_error(label + " is not set");
    }

    std::ifstream ifsVelOpt(velOptFile);
    if(!ifsVelOpt.is_open()) {
      throw std::runtime_error("Failed to open file: " + velOptFile);
    }
    std::vector<std::vector<double>> velOptTmp;

    while(getline(ifsVelOpt, str)) {
      std::istringstream iss(str);
      std::vector<double> tmp;
      for(int d = 0; d < dim; d++) {
        getline(iss, str, ' ');
        tmp.push_back(stod(str));
      }
      velOptTmp.push_back(tmp);
    }
    ifsVelOpt.close();
    velOpt[num] = velOptTmp;
    num++;
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readPostInverseFlowRateParameter()
{
  std::string str, base_label, label;
  base_label = "/PostInverseFrowRate";

  label = base_label + "/crossSection";
  if(!tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/flowRateVelDir";
  if(!tp.getInspectedValue(label, flowRateVelDir)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "xy") {
    crossSection = CrossSection::XY;
  } else if(str == "yz") {
    crossSection = CrossSection::YZ;
  } else if(str == "zx") {
    crossSection = CrossSection::ZX;
  }

  label = base_label + "/crossPoint";
  if(!tp.getInspectedValue(label, crossPoint)) {
    throw std::runtime_error(label + " is not set");
  }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readVoxelCreationParameter()
{
  std::string str, base_label, label;

  int tmpInt[dim];
  double tmpDouble[dim];

  base_label = "/VoxelCreation";

  label = base_label + "/inputDir";
  if(!tp.getInspectedValue(label, inputDir)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/timeMax";
  if(!tp.getInspectedValue(label, timeMax)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nSnapShot";
  if(!tp.getInspectedValue(label, nSnapShot)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/snapInterval";
  if(!tp.getInspectedValue(label, snapInterval)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/snapTimeBeginItr";
  if(!tp.getInspectedValue(label, snapTimeBeginItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nNodesInDataCell";
  if(!tp.getInspectedValue(label, nNodesInCellData)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/origin";
  if(!tp.getInspectedVector(label, tmpDouble, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  xOrigin = tmpDouble[0];
  yOrigin = tmpDouble[1];
  zOrigin = tmpDouble[2];

  label = base_label + "/nxData";
  if(!tp.getInspectedVector(label, tmpInt, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  nxData = tmpInt[0];
  nyData = tmpInt[1];
  nzData = tmpInt[2];

  label = base_label + "/lxData";
  if(!tp.getInspectedVector(label, tmpDouble, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  lxData = tmpDouble[0];
  lyData = tmpDouble[1];
  lzData = tmpDouble[2];

  dxData = lxData / (double)nxData;
  dyData = lyData / (double)nyData;
  dzData = lzData / (double)nzData;

  nDataCellsGlobal = nxData * nyData * nzData;

  label = base_label + "/nxOpt";
  if(!tp.getInspectedVector(label, tmpInt, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  nxOpt = tmpInt[0];
  nyOpt = tmpInt[1];
  nzOpt = tmpInt[2];

  label = base_label + "/lxOpt";
  if(!tp.getInspectedVector(label, tmpDouble, dim)) {
    throw std::runtime_error(label + " is not set");
  }

  lxOpt = tmpDouble[0];
  lyOpt = tmpDouble[1];
  lzOpt = tmpDouble[2];

  dxOpt = lxOpt / (double)nxOpt;
  dyOpt = lyOpt / (double)nyOpt;
  dzOpt = lzOpt / (double)nzOpt;

  nCellsOptGlobal = nxOpt * nxOpt * nzOpt;
  nNodesOptGlobal = (nxOpt + 1) * (nyOpt + 1) * (nzOpt + 1);
}
