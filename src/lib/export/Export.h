/**
 * @file Export.h
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef EXPORT_H
#define EXPORT_H

#include <iostream>
#include "Array.h"
#include "DataGrid.h"
#include "Node.h"
#include "Cell.h"
#include "MyMPI.h"

extern MyMPI mpi;

class EXPORT
{
public:
  // DAT
  template <typename T>
  static void exportScalarDataDAT(const std::string &file, const std::vector<T> &vec);
  template <typename T>
  static void exportVectorDataDAT(const std::string &file, const std::vector<std::vector<T>> &vec);
  template <typename T>
  static void exportMapDataDAT(const std::string &file, const std::map<int, std::vector<T>> &dataMap);
  template <typename T>
  static void exportMapDataDAT(const std::string &file, const std::map<int, T> &dataMap);

  static void exportCellDataDAT(const std::string &file, Cell &cell);
  static void exportNodeDataDAT(const std::string &file, Node &node);

  // BIN
  template <typename T>
  static void exportScalarDataBIN(const std::string &file, std::vector<T> &vec);
  template <typename T>
  static void exportVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec);

  // VTI
  // takes std::vector argument
  template <typename T>
  static void exportScalarPointDataVTI(const std::string &file, const char *dataName, std::vector<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  template <typename T>
  static void exportVectorPointDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<T>> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  template <typename T>
  static void exportScalarCellDataVTI(const std::string &file, const char *dataName, std::vector<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  template <typename T>
  static void exportVectorCellDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<T>> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);

  // takes custom array argument
  template <typename T>
  static void exportScalarPointDataVTI(const std::string &file, const char *dataName, Array1D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  template <typename T>
  static void exportVectorPointDataVTI(const std::string &file, const char *dataName, Array2D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  template <typename T>
  static void exportScalarCellDataVTI(const std::string &file, const char *dataName, Array1D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  template <typename T>
  static void exportVectorCellDataVTI(const std::string &file, const char *dataName, Array2D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);

  static void exportVelocityDataVTI(const std::string &file, DataGrid &data, const int t);

  // VTU
  // takes std::vector argument
  template <typename T>
  static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<T> &p);
  template <typename T>
  static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<T>> &p);
  template <typename T>
  static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<T> &c);
  template <typename T>
  static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<T>> &c);

  // takes custom array argument
  template <typename T>
  static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array1D<T> &p);
  template <typename T>
  static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<T> &p);
  template <typename T>
  static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array1D<T> &c);
  template <typename T>
  static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<T> &c);

  static void exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell);
  static void exportPhiVTU(const std::string &file, Node &node, Cell &cell);
};

#include "ExportDAT.inl"
#include "ExportBIN.inl"
#include "ExportVTU.inl"
#include "ExportVTI.inl"

#endif // EXPORT_H