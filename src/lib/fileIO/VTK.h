/**
 * @file VTK.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef VTK_H
#define VTK_H

#include <iostream>
#include "Array.h"
#include "Node.h"
#include "Cell.h"
#include "MyMPI.h"
#include "DataGrid.h"

extern MyMPI mpi;

class VTK
{
    public:
        static void exportScalarPointDataVTI(const std::string &file, const char *dataName, std::vector<double> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
        static void exportVectorPointDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<double>> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
        static void exportScalarCellDataVTI(const std::string &file, const char *dataName, std::vector<double> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
        static void exportVectorCellDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<double>> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
        static void exportVelocityDataVTI(const std::string &file, DataGrid &data, const int t);
     
        static void exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<double> &p);
        static void exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<double>> &p);
        static void exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<double> &c);
        static void exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<double>> &c);
        static void exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell);
        static void exportPhiVTU(const std::string &file, Node &node, Cell &cell);   
};

class BIN
{
    public:
        static void exportScalarDataBIN(const std::string &file, std::vector<double> &vec);
        static void exportVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec);
        static void importScalarDataBIN(const std::string &file, std::vector<double> &vec);
        static void importVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec);

};

#endif