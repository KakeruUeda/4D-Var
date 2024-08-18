/**
 * @file Export.h
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef EXPORT_H
#define EXPORT_H

#include "Array.h"
#include "Cell.h"
#include "MyMPI.h"
#include "Node.h"
#include <iostream>
#include <fstream>

extern MyMPI mpi;

class EXPORT
{
public:
  // DAT
  template <typename T>
  static void exportScalarDataDAT(const std::string &file, const std::vector<T> &vec);
  template <typename T>
  static void exportScalarDataDAT(const std::string &file, const std::set<T> &set);
  template <typename T>
  static void exportVectorDataDAT(const std::string &file, const std::vector<std::vector<T>> &vec);
  template <typename T>
  static void exportMapScalarDataDAT(const std::string &file, const std::map<int, T> &dataMap);
  template <typename T>
  static void exportMapVectorDataDAT(const std::string &file, const std::map<int, std::vector<T>> &dataMap);

  // BIN
  template <typename T>
  static void exportScalarDataBIN(const std::string &file, const std::vector<T> &vec);
  template <typename T>
  static void exportVectorDataBIN(const std::string &file, const std::vector<std::vector<T>> &vec);

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
  
  // and takes time argument
  template <typename T>
  static void exportScalarPointDataVTI(const std::string &file, const char *dataName, Array2D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz, const int t);
  template <typename T>
  static void exportVectorPointDataVTI(const std::string &file, const char *dataName, Array3D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz, const int t);
  template <typename T>
  static void exportScalarCellDataVTI(const std::string &file, const char *dataName, Array2D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz, const int t);
  template <typename T>
  static void exportVectorCellDataVTI(const std::string &file, const char *dataName, Array3D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz, const int t);

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
  
  // and takes time argument
  template <typename T>
  static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<T> &p, const int t);
  template <typename T>
  static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array3D<T> &p, const int t);
  template <typename T>
  static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<T> &c, const int t);
  template <typename T>
  static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array3D<T> &c, const int t);
  
  // These are ramdom export functions
  static void exportCellDataDAT(const std::string &file, Cell &cell);
  static void exportNodeDataDAT(const std::string &file, Node &node);
  static void exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell);
  static void exportPhiVTU(const std::string &file, Node &node, Cell &cell);
};

#include "ExportBIN.inl"
#include "ExportDAT.inl"
#include "ExportVTI.inl"
#include "ExportVTU.inl"

#endif  // EXPORT_H