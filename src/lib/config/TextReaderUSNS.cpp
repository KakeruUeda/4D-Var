/**
 * @file TextReaderUSNS.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderUSNS::readGridInfo()
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

  if(ptr->gridType == GridType::STRUCTURED) {
    sub_label = base_label + "/StructuredGrid";
    label = sub_label + "/nx";

    if(!ptr->tp.getInspectedVector(label, tmpInt, 3)) {
      throw std::runtime_error(label + " is not set.");
    }

    ptr->nx = tmpInt[0];
    ptr->ny = tmpInt[1];
    ptr->nz = tmpInt[2];

    ptr->nStrCellsGlobal = ptr->nx * ptr->ny * ptr->nz;
    ptr->nStrNodesGlobal = (ptr->nx + 1) * (ptr->ny + 1) * (ptr->nz + 1);

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

    label = sub_label + "/isOnlyFluidGrid";
    if(!ptr->tp.getInspectedValue(label, str)) {
      throw std::runtime_error(label + " is not set");
    }

    if(str == "YES") {
      ptr->isOnlyFluidGrid = true;
    } else if(str == "NO") {
      ptr->isOnlyFluidGrid = false;
    } else {
      throw std::runtime_error("Yes or No is not set");
    }

    if(ptr->isOnlyFluidGrid) {
      label = sub_label + "/fluidUniqueNodes";
      std::string fluidGridFile;
      if(!ptr->tp.getInspectedValue(label, fluidGridFile)) {
        throw std::runtime_error(label + " is not set");
      }
      IMPORT::importScalarDataDAT<int>(fluidGridFile, ptr->vecFluidUniqueNodes);
    }
  }

  sub_label = base_label + "/BaseGrid";
  label = sub_label + "/node";

  std::string nodeFile;
  if(!ptr->tp.getInspectedValue(label, nodeFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorDataDAT<double>(nodeFile, ptr->node);
  ptr->nNodesGlobal = ptr->node.size();

  std::string cellFile;
  label = sub_label + "/cell";
  if(!ptr->tp.getInspectedValue(label, cellFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorDataDAT<int>(cellFile, ptr->cell);
  ptr->nCellsGlobal = ptr->cell.size();

  std::string nodeIdFile;
  std::string cellIdFile;

  sub_label = base_label + "/SubGrid";

  label = sub_label + "/nodeId";
  if(!ptr->tp.getInspectedValue(label, nodeIdFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<int>(nodeIdFile, ptr->nodeId);

  label = sub_label + "/cellId";
  if(!ptr->tp.getInspectedValue(label, cellIdFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importScalarDataDAT<int>(cellIdFile, ptr->cellId);
}

/*****************************
 * @brief Read text parameter.
 */
void Config::TextReaderUSNS::readBoundaryInfo()
{
  std::string str, base_label, label;
  std::string nodeFile;

  base_label = "/Boundary";

  std::vector<int> vkey; 
  std::vector<std::vector<double>> vval;

  std::string velFile, preFile;

  label = base_label + "/velocityDirichlet";
  if(!ptr->tp.getInspectedValue(label, velFile)) {
    throw std::runtime_error(label + " is not set");
  }

  IMPORT::importVectorMapDataDAT<int, double>(velFile, vkey, vval);

  if(vkey.size() != vval.size()) {
    throw std::runtime_error("Size of vi and vvi is different");
  }

  for(int i = 0; i < vkey.size(); i++) {
    ptr->vDirichlet[vkey[i]] = vval[i];
  }

  label = base_label + "/pressureDirichlet";
  if(!ptr->tp.getInspectedValue(label, preFile)) {
    throw std::runtime_error(label + " is not set");
  }

  std::vector<int> pkey; 
  std::vector<double> pval;
  IMPORT::importScalarMapDataDAT<int, double>(preFile, pkey, pval);

  if(pkey.size() != pval.size()) {
    throw std::runtime_error("Size of vi and vd is different");
  }

  for(int i = 0; i < pkey.size(); i++) {
    ptr->pDirichlet[pkey[i]] = pval[i];
  }
}