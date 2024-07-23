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

enum class DataType
{
    VELOCITY = 0,
    PRESSURE = 1,
    FEEDBACK = 2,
    ADJOINT_W = 3,
    ADJOINT_Q = 4,
    ADJOINT_L = 5
};

class VTK
{
    public:
        // VTI
        void exportDataVTI(const std::string file, DataGrid &data, const int &t, const int &dim);
        void exportVelocityDataVTI(const std::string file, DataGrid &data, const int &t, const int &dim);
        void exportSolutionVTI(const string &file, std::vector<std::vector<double>> &v, std::vector<double> &p,
                               const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
        void exportAdjointSolutionVTI(const string &file, std::vector<std::vector<double>> &w, std::vector<double> &q, std::vector<std::vector<double>> &l,
                                      const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
        void exportNodeVTI(const string &file, std::vector<std::vector<double>> &node,
                           const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);

        // VTU
        void exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell);
        void exportPhiVTU(const std::string &file, Node &node, Cell &cell);
        void exportSolutionVTU(const std::string &file, Node &node, Cell &cell, DataType dataType);
        void exportMainVariablesVTU(const std::string &file, Node &node, Cell &cell, const int t, DataType dataType);
        void exportSnapShotVTU(const std::string &file, Node &node, Cell &cell, SnapShot &snap, const int &snapCount);
        void exportAdjointSolutionVTU(const std::string &file, Node &node, Cell &cell, DataType dataType);
        void exportFeedbackForceVTU(const std::string &file, Node &node, Cell &cell, std::vector<std::vector<double>> &feedback);
};


#endif