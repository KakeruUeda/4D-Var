/**
 * @file FileIO.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef FILEIO_H
#define FILEIO_H

#include <iostream>
#include "Array.h"
#include "DataGrid.h"
#include "Node.h"
#include "Cell.h"
#include "MyMPI.h"

extern MyMPI mpi;

class DAT
{
public:
  template <typename T>
  static void exportScalarDataDAT(const std::string &file, const std::vector<T> &vec);
  template <typename T>
  static void exportVectorDataDAT(const std::string &file, const std::vector<std::vector<T>> &vec);
  template <typename T>
  static void importScalarDataDAT(const std::string &file, std::vector<T> &vec);
  template <typename T>
  static void importVectorDataDAT(const std::string &file, std::vector<std::vector<T>> &vec);
  template <typename T>
  static void exportMapDataDAT(const std::string &file, const std::map<int, std::vector<T>> &dataMap);
  template <typename T>
  static void exportMapDataDAT(const std::string &file, const std::map<int, T> &dataMap);

  static void exportCellDataDAT(const std::string &file, Cell &cell);
  static void exportNodeDataDAT(const std::string &file, Node &node);
};

class BIN
{
public:
  static void exportScalarDataBIN(const std::string &file, std::vector<double> &vec);
  static void exportVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec);
  static void importScalarDataBIN(const std::string &file, std::vector<double> &vec);
  static void importVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec);
};

class VTK
{
public:
  // takes std::vector argument (vti)
  static void exportScalarPointDataVTI(const std::string &file, const char *dataName, std::vector<double> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  static void exportVectorPointDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<double>> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  static void exportScalarCellDataVTI(const std::string &file, const char *dataName, std::vector<double> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  static void exportVectorCellDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<double>> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  
  // takes custom array argument (vti)
  static void exportScalarPointDataVTI(const std::string &file, const char *dataName, Array1D<double> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  static void exportVectorPointDataVTI(const std::string &file, const char *dataName, Array2D<double> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  static void exportScalarCellDataVTI(const std::string &file, const char *dataName, Array1D<double> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
  static void exportVectorCellDataVTI(const std::string &file, const char *dataName, Array2D<double> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);

  static void exportVelocityDataVTI(const std::string &file, DataGrid &data, const int t);

  // takes std::vector argument (vtu)
  static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<double> &p);
  static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<double>> &p);
  static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<double> &c);
  static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<double>> &c);
  
  // takes custom array argument (vtu)
  static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array1D<double> &p);
  static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<double> &p);
  static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array1D<double> &c);
  static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<double> &c);

  static void exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell);
  static void exportPhiVTU(const std::string &file, Node &node, Cell &cell);
};

class BINTMP
{
public:
    template <typename T>
    static void exportScalarDataBIN(const std::string &file, std::vector<T> &vec);
    template <typename T>
    static void exportVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec);
    template <typename T>
    static void importScalarDataBIN(const std::string &file, std::vector<T> &vec);
    template <typename T>
    static void importVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec);
};

class VTKTMP
{
public:
    // takes std::vector argument (vti)
    template <typename T>
    static void exportScalarPointDataVTI(const std::string &file, const char *dataName, std::vector<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    template <typename T>
    static void exportVectorPointDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<T>> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    template <typename T>
    static void exportScalarCellDataVTI(const std::string &file, const char *dataName, std::vector<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    template <typename T>
    static void exportVectorCellDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<T>> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);

    // takes custom array argument (vti)
    template <typename T>
    static void exportScalarPointDataVTI(const std::string &file, const char *dataName, Array1D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    template <typename T>
    static void exportVectorPointDataVTI(const std::string &file, const char *dataName, Array2D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    template <typename T>
    static void exportScalarCellDataVTI(const std::string &file, const char *dataName, Array1D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    template <typename T>
    static void exportVectorCellDataVTI(const std::string &file, const char *dataName, Array2D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);

    static void exportVelocityDataVTI(const std::string &file, DataGrid &data, const int t);

    // takes std::vector argument (vtu)
    template <typename T>
    static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<T> &p);
    template <typename T>
    static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<T>> &p);
    template <typename T>
    static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<T> &c);
    template <typename T>
    static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<T>> &c);
    
    // takes custom array argument (vtu)
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

#include "DAT.inl"
#include "BIN.inl"
#include "VTU.inl"
#include "VTI.inl"

#endif