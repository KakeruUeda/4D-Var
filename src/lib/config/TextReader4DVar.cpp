/**
 * @file TextReader4DVar.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/**
 * @brief Read text parameter.
 */
void TextReader4DVar::readBasicInfo(Config &conf)
{
  TextReaderInterface::readBasicInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReader4DVar::readBoundaryInfo(Config &conf)
{
  TextReaderInterface::readBoundaryInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReader4DVar::readPhysicalInfo(Config &conf)
{
  TextReaderInterface::readPhysicalInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReader4DVar::readDarcyInfo(Config &conf)
{
  TextReaderInterface::readDarcyInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReader4DVar::readTimeInfo(Config &conf)
{
  TextReaderInterface::readTimeInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReader4DVar::readGridInfo(Config &conf)
{
  TextReaderInterface::readGridInfo(conf);

  std::string str, base_label, sub_label, label;

  base_label = "/Grid";

  if(conf.gridType == GridType::STRUCTURED) {
    conf.nStrCellsGlobal = conf.nx * conf.ny * conf.nz;
    conf.nStrNodesGlobal = (conf.nx + 1) * (conf.ny + 1) * (conf.nz + 1);

    sub_label = base_label + "/StructuredGrid";
    label = sub_label + "/image";

    std::string imageFile;
    if(!conf.tp.getInspectedValue(label, imageFile)) {
      throw std::runtime_error(label + " is not set");
    }

    IMPORT::importScalarDataDAT<double>(imageFile, conf.phi);

    label = sub_label + "/fluidExtraction";
    if(!conf.tp.getInspectedValue(label, str)) {
      throw std::runtime_error(label + " is not set");
    }

    if(str == "ON") {
      conf.fluidExtraction = ON;
    } else if(str == "OFF") {
      conf.fluidExtraction = OFF;
    } else {
      throw std::runtime_error("ON or OFF is not set");
    }

    if(conf.fluidExtraction) {
      label = sub_label + "/fluidUniqueNodes";
      std::string fluidGridFile;
      if(!conf.tp.getInspectedValue(label, fluidGridFile)) {
        throw std::runtime_error(label + " is not set");
      }
      IMPORT::importScalarDataDAT<int>(fluidGridFile, conf.vecFluidUniqueNodes);
    } else {
      conf.vecFluidUniqueNodes.resize(0);
    }
  }

  sub_label = base_label + "/BaseGrid";
  label = sub_label + "/node";

  std::string nodeFile;
  if(!conf.tp.getInspectedValue(label, nodeFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorDataDAT<double>(nodeFile, conf.node);
  conf.nNodesGlobal = conf.node.size();

  std::string cellFile;
  label = sub_label + "/cell";
  if(!conf.tp.getInspectedValue(label, cellFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorDataDAT<int>(cellFile, conf.cell);
  conf.nCellsGlobal = conf.cell.size();

  std::string nodeIdFile;
  std::string cellIdFile;

  sub_label = base_label + "/SubGrid";

  label = sub_label + "/nodeId";
  if(!conf.tp.getInspectedValue(label, nodeIdFile)) {
    throw std::runtime_error(label + " is not set");
  }
  IMPORT::importScalarDataDAT<int>(nodeIdFile, conf.nodeId);

  label = sub_label + "/cellId";
  if(!conf.tp.getInspectedValue(label, cellIdFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<int>(cellIdFile, conf.cellId);
}



void TextReader4DVar::readInverseInfo(Config &conf)
{
  std::string str, base_label, label;
  base_label = "/Inverse";

  label = base_label + "/aCF";
  if(!conf.tp.getInspectedValue(label, conf.aCF)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/bCF";
  if(!conf.tp.getInspectedValue(label, conf.bCF)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/gCF";
  if(!conf.tp.getInspectedValue(label, conf.gCF)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/alphaX0";
  if(!conf.tp.getInspectedValue(label, conf.alphaX0)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/alphaX";
  if(!conf.tp.getInspectedValue(label, conf.alphaX)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/outputItr";
  if(!conf.tp.getInspectedValue(label, conf.outputItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/loopMax";
  if(!conf.tp.getInspectedValue(label, conf.loopMax)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/controlBoundary";

  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.planeDir.resize(2, 0);

  if(str == "left") {
    conf.controlBoundary = ControlBoundaryFace::left;
    conf.planeDir[0] = 1, conf.planeDir[1] = 2;
  } else if(str == "right") {
    conf.controlBoundary = ControlBoundaryFace::right;
    conf.planeDir[0] = 1, conf.planeDir[1] = 2;
  } else if(str == "top") {
    conf.controlBoundary = ControlBoundaryFace::top;
    conf.planeDir[0] = 0, conf.planeDir[1] = 2;
  } else if(str == "bottom") {
    conf.controlBoundary = ControlBoundaryFace::bottom;
    conf.planeDir[0] = 0, conf.planeDir[1] = 2;
  } else if(str == "front") {
    conf.controlBoundary = ControlBoundaryFace::front;
    conf.planeDir[0] = 0, conf.planeDir[1] = 1;
  } else if(str == "back") {
    conf.controlBoundary = ControlBoundaryFace::back;
    conf.planeDir[0] = 0, conf.planeDir[1] = 1;
  }

  label = base_label + "/nControlNodesInCell";
  if(!conf.tp.getInspectedValue(label, conf.nControlNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }

  std::string controlBoundaryNodeMap;
  label = base_label + "/controlBoundaryNodeMap";

  if(!conf.tp.getInspectedValue(label, controlBoundaryNodeMap)) {
    throw std::runtime_error(label + " is not set");
  }
  IMPORT::importScalarDataDAT<int>(controlBoundaryNodeMap, conf.CBNodeMap);

  std::string controlBoundaryCellMap;
  label = base_label + "/controlBoundaryCellMap";

  if(!conf.tp.getInspectedValue(label, controlBoundaryCellMap)) {
    throw std::runtime_error(label + " is not set");
  }
  IMPORT::importScalarDataDAT<int>(controlBoundaryCellMap, conf.CBCellMap);

  std::string controlBoundaryNodeMapInCell;
  label = base_label + "/controlBoundaryNodeMapInCell";

  if(!conf.tp.getInspectedValue(label, controlBoundaryNodeMapInCell)) {
    throw std::runtime_error(label + " is not set");
  }
  IMPORT::importVectorDataDAT<int>(controlBoundaryNodeMapInCell, conf.CBNodeMapInCell);
}

void TextReader4DVar::readDataInfo(Config &conf)
{
  std::string str, base_label, label;
  base_label = "/Data";

  int tmpInt[3] = {0, 0, 0};
  double tmpDouble[3] = {0.0, 0.0, 0.0};

  label = base_label + "/nSnapShot";
  if(!conf.tp.getInspectedValue(label, conf.nSnapShot)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/snapInterval";
  if(!conf.tp.getInspectedValue(label, conf.snapInterval)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/snapTimeBeginItr";
  if(!conf.tp.getInspectedValue(label, conf.snapTimeBeginItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nNodesInDataCell";
  if(!conf.tp.getInspectedValue(label, conf.nNodesInCellData)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/nxData";
  if(!conf.tp.getInspectedVector(label, tmpInt, 3)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.nxData = tmpInt[0];
  conf.nyData = tmpInt[1];
  conf.nzData = tmpInt[2];

  label = base_label + "/lxData";
  if(!conf.tp.getInspectedVector(label, tmpDouble, 3)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.lxData = tmpDouble[0];
  conf.lyData = tmpDouble[1];
  conf.lzData = tmpDouble[2];

  conf.dxData = conf.lxData / (double)conf.nxData;
  conf.dyData = conf.lyData / (double)conf.nyData;
  conf.dzData = conf.lzData / (double)conf.nzData;

  conf.nDataCellsGlobal = conf.nxData * conf.nyData * conf.nzData;

  label = base_label + "/voxelVelocity";

  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "average") {
    conf.vvox = VoxelVelocity::AVERAGE;
  } else if(str == "interpolation") {
    conf.vvox = VoxelVelocity::INTERPOLATION;
  } else {
    throw std::runtime_error("undefined voxelVelocity");
  }

  label = base_label + "/dataDir";
  if(!conf.tp.getInspectedValue(label, conf.dataDir)) {
    throw std::runtime_error(label + " is not set");
  }
}
