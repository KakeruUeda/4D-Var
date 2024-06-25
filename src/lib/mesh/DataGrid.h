#ifndef DATAGRID_H
#define DATAGRID_H

#include <iostream>
#include "metis.h"
#include "Config.h"
#include "Cell.h"
#include "Node.h"
#include "Config.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include "MathFEM.h"

struct VoxelInfo
{
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> vCFD;
    std::vector<std::vector<double>> vMRI;
    std::vector<double> center;
    std::vector<int> cellChildren;

    void setNearCell(Node &node, Cell &cell, const double &length, const int &dim);
    void averageVelocity(Cell &cell, std::vector<std::vector<double>> &_v, 
                         const int &t, const int &nNodesInCell, const int &dim);
    void gaussIntegral(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent, 
                       std::vector<std::vector<double>> &velCurrent, double &weightIntegral, 
                       const int &nNodesInCell, const double &detJ, const double &weight, 
                       const int &t, const int &dim);
};

class DataGrid
{
    public:
        DataGrid(){}  
        DataGrid(Config &conf);

        int nx, ny, nz;
        double dx, dy, dz; 
        double lx, ly, lz;
        double xOrigin, yOrigin, zOrigin;

        int nCellsGlobal;
        int nNodesGlobal;
        int nNodesInCell;

        int nSnapShot;
        int snapInterval;
        
        inline VoxelInfo& operator()(int x)
        { return data[x]; }
        inline VoxelInfo& operator()(int y, int x)
        { return data[y * nx + x]; }
        inline VoxelInfo& operator()(int z, int y, int x)
        { return data[z * nx * ny + y * nx + x]; }

        inline int size()
        { return data.size(); }
        inline void resize(int n)
        { data.resize(n); }

        std::vector<VoxelInfo> data;
};


#endif