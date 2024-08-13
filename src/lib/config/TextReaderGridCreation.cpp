/**
 * @file Config.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void Config::TextReaderGridCreation::readGridInfo()
{
  std::string str, base_label, sub_label, label;
  int tmpInt[3];
  double tmpDouble[3];

  base_label = "/Grid";
  sub_label = base_label + "/GridType";

  label = sub_label + "/type";
  if(!ptr->tp.getInspectedValue(label, str)) {
    throw std::runtime_error(label + " is not set");
  }

  if(str == "Structured") {
    ptr->gridType = GridType::STRUCTURED;
  } else if(str == "Unstructured") {
    ptr->gridType = GridType::UNSTRUCTURED;
  } else {
    throw std::runtime_error("Unknown GridType");
  }

  label = sub_label + "/nNodesInCell";
  if(!ptr->tp.getInspectedValue(label, ptr->nNodesInCell)) {
    throw std::runtime_error(label + " is not set");
  }

  sub_label = base_label + "/StructuredGrid";
  label = sub_label + "/nx";

  if(!ptr->tp.getInspectedVector(label, tmpInt, 3)) {
    throw std::runtime_error(label + " is not set.");
  }

  ptr->nx = tmpInt[0];
  ptr->ny = tmpInt[1];
  ptr->nz = tmpInt[2];

  ptr->nCellsGlobal = ptr->nx * ptr->ny * ptr->nz;
  ptr->nNodesGlobal = (ptr->nx + 1) * (ptr->ny + 1) * (ptr->nz + 1);

  label = sub_label + "/lx";
  if(!ptr->tp.getInspectedVector(label, tmpDouble, ptr->dim)) {
    throw std::runtime_error(label + " is not set");
  }

  ptr->lx = tmpDouble[0];
  ptr->ly = tmpDouble[1];
  ptr->lz = tmpDouble[2];

  ptr->dx = ptr->lx / (double)ptr->nx;
  ptr->dy = ptr->ly / (double)ptr->ny;
  ptr->dz = ptr->lz / (double)ptr->nz;

  label = sub_label + "/image";

  std::string imageFile;
  if(!ptr->tp.getInspectedValue(label, imageFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<double>(imageFile, ptr->phi);
  ptr->setStrGrid();  // define str grid
}

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderGridCreation::readStructuredBoundaryInfo()
{
  std::string base_label = "/StructuredBoundary";

  auto readBoundary = [&](const std::string &face)
  {
    std::string labelFace = base_label + "/" + face;
    std::string labelType = labelFace + "/type";
    std::string labelValue = labelFace + "/value";

    std::string bdTypeTmp;
    if(!ptr->tp.getInspectedValue(labelType, bdTypeTmp)) {
      throw std::runtime_error(labelType + " is not set");
    }

    if(bdTypeTmp == "v") {
      double value[3];
      if(!ptr->tp.getInspectedVector(labelValue, value, 3)) {
        throw std::runtime_error(labelValue + " is not set");
      }
      ptr->setBoundaryVelocityValue(face, value);
    } else if(bdTypeTmp == "p") {
      double value;
      if(!ptr->tp.getInspectedValue(labelValue, value)) {
        throw std::runtime_error(labelValue + " is not set");
      }
      ptr->setBoundaryPressureValue(face, value);
    } else if(bdTypeTmp == "poiseuille") {
      std::string label = labelFace + "/center";
      if(!ptr->tp.getInspectedVector(label, ptr->center, 3)) {
        throw std::runtime_error(label + " is not set");
      }

      label = labelFace + "/R";
      if(!ptr->tp.getInspectedValue(label, ptr->R)) {
        throw std::runtime_error(label + " is not set");
      }

      label = labelFace + "/Q";
      if(!ptr->tp.getInspectedValue(label, ptr->Q)) {
        throw std::runtime_error(label + " is not set");
      }
      ptr->setBoundaryPoiseuilleValue(face);
    } else if(bdTypeTmp == "free") {
      // "free" type specific logic here, if any
    } else {
      throw std::runtime_error("label " + bdTypeTmp + " undefined");
    }
  };

  std::vector<std::string> faces = {"bottom", "top", "left", "right", "front", "back"};
  for(const auto &face : faces) {
    readBoundary(face);
  }
}
