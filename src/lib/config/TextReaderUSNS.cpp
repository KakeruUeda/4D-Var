/**
 * @file TextReaderUSNS.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/*****************************
 * @brief Read text parameter.
 */
void TextReaderUSNS::readBasicInfo(Config &conf)
{
  TextReaderInterface::readBasicInfo(conf);
}

/*****************************
 * @brief Read text parameter.
 */
void TextReaderUSNS::readBoundaryInfo(Config &conf)
{
  TextReaderInterface::readBoundaryInfo(conf);
}

/*****************************
 * @brief Read text parameter.
 */
void TextReaderUSNS::readPhysicalInfo(Config &conf)
{
  TextReaderInterface::readPhysicalInfo(conf);
}

/*****************************
 * @brief Read text parameter.
 */
void TextReaderUSNS::readDarcyInfo(Config &conf)
{
  TextReaderInterface::readDarcyInfo(conf);
}

/*****************************
 * @brief Read text parameter.
 */
void TextReaderUSNS::readTimeInfo(Config &conf)
{
  TextReaderInterface::readTimeInfo(conf);
}

/*****************************
 * @brief Read text parameter.
 */
void TextReaderUSNS::readGridInfo(Config &conf)
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
    }else{
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

