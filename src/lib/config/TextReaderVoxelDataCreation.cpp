/**
 * @file TextReaderUSNS.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/**
 * @brief Read text parameter.
 */
void TextReaderVoxelDataCreation::readBasicInfo(Config &conf)
{
  TextReaderInterface::readBasicInfo(conf);
}

/**
 * @brief Read text parameter.
 */
void TextReaderVoxelDataCreation::readSnapInfo(Config &conf)
{
  std::string str, base_label, label;

  base_label = "/SnapShot";

  label = base_label + "/timeMax";
  if(!conf.tp.getInspectedValue(label, conf.timeMax)) {
    throw std::runtime_error(label + " is not set");
  }

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
}

/**
 * @brief Read text parameter.
 */
void TextReaderVoxelDataCreation::readGridInfo(Config &conf)
{
  TextReaderInterface::readGridInfo(conf);

  std::string str, base_label, sub_label, label;
  double tmpDouble[3];
  int tmpInt[3];

  base_label = "/Grid";

  if(conf.gridType == GridType::STRUCTURED) {
    conf.nCellsGlobal = conf.nx * conf.ny * conf.nz;
    conf.nNodesGlobal = (conf.nx + 1) * (conf.ny + 1) * (conf.nz + 1);
    conf.vecFluidUniqueNodes.resize(0);
    conf.setStrGrid();
  } else if(conf.gridType == GridType::UNSTRUCTURED) {
    throw std::runtime_error("Unstructured grid is not supported yet");
  }

  sub_label = base_label + "/DataGrid";

  label = sub_label + "/numerical_velocity_space";

  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "weighted_average") {
    conf.vel_space = VoxelVelocity::WEIGHTED_AVERAGE;
  }else if(str == "average") {
    conf.vel_space = VoxelVelocity::AVERAGE;
  } else if(str == "interpolation") {
    conf.vel_space = VoxelVelocity::INTERPOLATION;
  } else {
    throw std::runtime_error("undefined numerical_velocity_space");
  }

  label = sub_label + "/numerical_velocity_time";

  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "weighted_average") {
    conf.vel_time = VoxelVelocity::WEIGHTED_AVERAGE;
  }else if(str == "average") {
    conf.vel_time = VoxelVelocity::AVERAGE;
  } else if(str == "interpolation") {
    conf.vel_time = VoxelVelocity::INTERPOLATION;
  } else {
    throw std::runtime_error("undefined numerical_velocity_time");
  }

  label = sub_label + "/nNodesInDataCell";
  if(!conf.tp.getInspectedValue(label, conf.nNodesInCellData)) {
    throw std::runtime_error(label + " is not set");
  }

  label = sub_label + "/nxData";
  if(!conf.tp.getInspectedVector(label, tmpInt, conf.dim)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.nxData = tmpInt[0];
  conf.nyData = tmpInt[1];
  conf.nzData = tmpInt[2];

  label = sub_label + "/lxData";
  if(!conf.tp.getInspectedVector(label, tmpDouble, conf.dim)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.lxData = tmpDouble[0];
  conf.lyData = tmpDouble[1];
  conf.lzData = tmpDouble[2];

  conf.dxData = conf.lxData / (double)conf.nxData;
  conf.dyData = conf.lyData / (double)conf.nyData;
  conf.dzData = conf.lzData / (double)conf.nzData;

  conf.nDataCellsGlobal = conf.nxData * conf.nyData * conf.nzData;

  sub_label = base_label + "/OptGrid";

  label = sub_label + "/origin";
  if(!conf.tp.getInspectedVector(label, tmpDouble, conf.dim)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.xOrigin = tmpDouble[0];
  conf.yOrigin = tmpDouble[1];
  conf.zOrigin = tmpDouble[2];

  label = sub_label + "/nxOpt";
  if(!conf.tp.getInspectedVector(label, tmpInt, conf.dim)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.nxOpt = tmpInt[0];
  conf.nyOpt = tmpInt[1];
  conf.nzOpt = tmpInt[2];

  label = sub_label + "/lxOpt";
  if(!conf.tp.getInspectedVector(label, tmpDouble, conf.dim)) {
    throw std::runtime_error(label + " is not set");
  }

  conf.lxOpt = tmpDouble[0];
  conf.lyOpt = tmpDouble[1];
  conf.lzOpt = tmpDouble[2];

  conf.dxOpt = conf.lxOpt / (double)conf.nxOpt;
  conf.dyOpt = conf.lyOpt / (double)conf.nyOpt;
  conf.dzOpt = conf.lzOpt / (double)conf.nzOpt;

  conf.nCellsOptGlobal = conf.nxOpt * conf.nxOpt * conf.nzOpt;
  conf.nNodesOptGlobal = (conf.nxOpt + 1) * (conf.nyOpt + 1) * (conf.nzOpt + 1);
}

/**
 * @brief Read text parameter.
 */
void TextReaderVoxelDataCreation::readOriginalInfo(Config &conf)
{
  std::string str, base_label, sub_label, label;

  base_label = "/Original";

  label = base_label + "/inputDir";
  if(!conf.tp.getInspectedValue(label, conf.inputDir)) {
    throw std::runtime_error(label + " is not set");
  }
}
