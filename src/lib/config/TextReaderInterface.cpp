/**
 * @file Config.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/**
 * @brief Read text parameter.
 */
void TextReaderInterface::readBasicInfo(Config &conf)
{
  std::string str, base_label, label;
  base_label = "/Base";
  label = base_label + "/dim";

  if(!conf.tp.getInspectedValue(label, conf.dim)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/outputDir";
  if(!conf.tp.getInspectedValue(label, conf.outputDir)) {
    throw std::runtime_error(label + " is not set");
  }
}

/** 
 * @brief Read interface grid info.
 */
void TextReaderInterface::readGridInfo(Config &conf)
{
  std::string str, base_label, sub_label, label;
  int tmpInt[3];
  double tmpDouble[3];

  base_label = "/Grid";
  sub_label = base_label + "/GridType";

  label = sub_label + "/type";
  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "Structured") {
    conf.gridType = GridType::STRUCTURED;
  } else if(str == "Unstructured") {
    conf.gridType = GridType::UNSTRUCTURED;
  } else {
    throw std::runtime_error("Unknown GridType");
  }

  label = sub_label + "/nNodesInCell";
  if(!conf.tp.getInspectedValue(label, conf.nNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }

  if(conf.gridType == GridType::STRUCTURED) {
    sub_label = base_label + "/StructuredGrid";
    label = sub_label + "/nx";

    if(!conf.tp.getInspectedVector(label, tmpInt, 3)) {
      throw std::runtime_error(label + " is not set.");
    }

    conf.nx = tmpInt[0];
    conf.ny = tmpInt[1];
    conf.nz = tmpInt[2];

    label = sub_label + "/lx";
    if(!conf.tp.getInspectedVector(label, tmpDouble, conf.dim)) {
      throw std::runtime_error(label + " is not set");
    }

    conf.lx = tmpDouble[0];
    conf.ly = tmpDouble[1];
    conf.lz = tmpDouble[2];

    conf.dx = conf.lx / (double)conf.nx;
    conf.dy = conf.ly / (double)conf.ny;
    conf.dz = conf.lz / (double)conf.nz;

  } else if(conf.gridType == GridType::UNSTRUCTURED) {
  } else {
    throw std::runtime_error("Unknown GridType");
  }
}

/**
 * @brief Read text parameter.
 */
void TextReaderInterface::readBoundaryInfo(Config &conf)
{
  std::string str, base_label, label;
  std::string nodeFile;

  base_label = "/Boundary";

  std::vector<int> vkey;
  std::vector<std::vector<double>> vval;

  std::string velFile, preFile;

  label = base_label + "/velocityDirichlet";
  if(!conf.tp.getInspectedValue(label, velFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorMapDataDAT<int, double>(velFile, vkey, vval);

  if(vkey.size() != vval.size()) {
    throw std::runtime_error("Size of vi and vvi is different");
  }

  for(int i = 0; i < vkey.size(); i++) {
    conf.vDirichlet[vkey[i]] = vval[i];
  }

  label = base_label + "/pressureDirichlet";
  if(!conf.tp.getInspectedValue(label, preFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::vector<int> pkey;
  std::vector<double> pval;
  IMPORT::importScalarMapDataDAT<int, double>(preFile, pkey, pval);

  if(pkey.size() != pval.size()) {
    throw std::runtime_error("Size of vi and vd is different");
  }

  for(int i = 0; i < pkey.size(); i++) {
    conf.pDirichlet[pkey[i]] = pval[i];
  }
}

/**
 * @brief Read text parameter.
 */
void TextReaderInterface::readPhysicalInfo(Config &conf)
{
  std::string str, base_label, label;

  base_label = "/PysicalParameter";
  label = base_label + "/rho";
  if(!conf.tp.getInspectedValue(label, conf.rho)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/mu";
  if(!conf.tp.getInspectedValue(label, conf.mu)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/L";
  if(!conf.tp.getInspectedValue(label, conf.L)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.nu = conf.mu / conf.rho;
}

/**
 * @brief Read text parameter.
 */
void TextReaderInterface::readTimeInfo(Config &conf)
{
  std::string str, base_label, label;

  base_label = "/TimeParameter";
  label = base_label + "/dt";
  if(!conf.tp.getInspectedValue(label, conf.dt)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/timeMax";
  if(!conf.tp.getInspectedValue(label, conf.timeMax)) {
    throw std::runtime_error(label + " is not set");
  }

  std::string ON_OFF;

  label = base_label + "/pulsatileFlow";
  if(!conf.tp.getInspectedValue(label, ON_OFF)) {
    throw std::runtime_error(label + " is not set");
  }

  if(ON_OFF == "ON") {
    conf.pulsatileFlow = ON;
  } else if(ON_OFF == "OFF") {
    conf.pulsatileFlow = OFF;
  } else {
    throw std::runtime_error("ON or OFF is not set");
  }

  label = base_label + "/pulseBeginItr";
  if(!conf.tp.getInspectedValue(label, conf.pulseBeginItr)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/T";
  if(!conf.tp.getInspectedValue(label, conf.T)) {
    throw std::runtime_error(label + " is not set");
  }
}

/**
 * @brief Read text parameter.
 */
void TextReaderInterface::readDarcyInfo(Config &conf)
{
  std::string str, base_label, label;
  base_label = "/DarcyParameter";

  label = base_label + "/alpha";
  if(!conf.tp.getInspectedValue(label, conf.alpha)) {
    throw std::runtime_error(label + " is not set");
  }

  label = base_label + "/resistance";
  if(!conf.tp.getInspectedValue(label, conf.resistance)) {
    throw std::runtime_error(label + " is not set");
  }
}
