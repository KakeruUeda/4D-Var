/**
 * @file Config.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void TextReaderGridCreation::readBasicInfo(Config &conf)
{
  TextReaderInterface::readBasicInfo(conf);
}

void TextReaderGridCreation::readGridInfo(Config &conf)
{
  TextReaderInterface::readGridInfo(conf);
  conf.nCellsGlobal = conf.nx * conf.ny * conf.nz;
  conf.nNodesGlobal = (conf.nx + 1) * (conf.ny + 1) * (conf.nz + 1);

  std::string str, base_label, sub_label, label;
  int tmpInt[3];
  double tmpDouble[3];

  base_label = "/Grid";

  if(conf.gridType == GridType::STRUCTURED) {
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

    conf.setStrGrid();  // define str grid
  }
}

/*****************************
 * @brief Read text parameter.
 */
void TextReaderGridCreation::readStructuredBoundaryInfo(Config &conf)
{
  std::string str, base_label, label;
  base_label = "/StructuredBoundary";

  auto readBoundary = [&](const std::string &face) {
    std::string labelFace = base_label + "/" + face;
    std::string labelType = labelFace + "/type";
    std::string labelValue = labelFace + "/value";

    std::string bdTypeTmp;
    if(!conf.tp.getInspectedValue(labelType, bdTypeTmp)) {
      throw std::runtime_error(labelType + " is not set");
    }

    if(bdTypeTmp == "v") {
      double value[3];
      if(!conf.tp.getInspectedVector(labelValue, value, 3)) {
        throw std::runtime_error(labelValue + " is not set");
      }
      conf.setBoundaryVelocityValue(face, value);
    } else if(bdTypeTmp == "p") {
      double value;
      if(!conf.tp.getInspectedValue(labelValue, value)) {
        throw std::runtime_error(labelValue + " is not set");
      }
      conf.setBoundaryPressureValue(face, value);
    } else if(bdTypeTmp == "poiseuille") {
      std::string label = labelFace + "/center";
      if(!conf.tp.getInspectedVector(label, conf.center, 3)) {
        throw std::runtime_error(label + " is not set");
      }

      label = labelFace + "/R";
      if(!conf.tp.getInspectedValue(label, conf.R)) {
        throw std::runtime_error(label + " is not set");
      }

      label = labelFace + "/Q";
      if(!conf.tp.getInspectedValue(label, conf.Q)) {
        throw std::runtime_error(label + " is not set");
      }
      conf.setBoundaryPoiseuilleValue(face);
    } else if(bdTypeTmp == "free") {
    } else {
      throw std::runtime_error("label " + bdTypeTmp + " undefined");
    }
  };

  std::vector<std::string> faces = {"bottom", "top", "left", "right", "front", "back"};
  for(const auto &face : faces) {
    readBoundary(face);
  }

  label = base_label + "/controlBoundaryExtraction";

  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "ON") {
    conf.CBExtraction = ON;
  } else if(str == "OFF") {
    conf.CBExtraction = OFF;
    return;
  } else {
    throw std::runtime_error("ON or OFF is not set");
  }

  label = base_label + "/inletControlBoundary";

  if(!conf.tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "left") {
    conf.inletCB = ControlBoundary::left;
  }
  if(str == "right") {
    conf.inletCB = ControlBoundary::right;
  }
  if(str == "top") {
    conf.inletCB = ControlBoundary::top;
  }
  if(str == "bottom") {
    conf.inletCB = ControlBoundary::bottom;
  }
  if(str == "front") {
    conf.inletCB = ControlBoundary::front;
  }
  if(str == "back") {
    conf.inletCB = ControlBoundary::back;
  }

  conf.setControlBoundary();
}
