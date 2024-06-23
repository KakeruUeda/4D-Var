#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "metis.h"
#include "Config.h"
#include "Cell.h"
#include "Node.h"
#include "Array.h"
#include "Boundary.h"
#include "Output.h"
#include "Config.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include "MathFEM.h"
 
class Grid
{
    public:
        Grid(){}
        Grid(Config &conf);
        virtual ~Grid(){}

        GridType gridType;
        Cell cell; 
        Node node;
        DirichletBoundary dirichlet;
        OutputVTU outputVTU;
        
        int nx, ny, nz;
        double lx, ly, lz, dx, dy, dz;

        int dim;
        int extractFluid;
        int rowStart, rowEnd;

        int nNodesGlobal, nCellsGlobal, nDofsGlobal;
        int nNodesLocal, nCellsLocal, nDofsLocal;

        void setStructuredGrid(
            const int nxCells, const int nyCells, const int nzCells, 
            const int nxNodes, const int nyNodes, const int nzNodes, 
            const double dx, const double dy, const double dz,
            const int nNodesInCell, const int dim, 
            Cell &cell, Node &node);
        void divideWholeGrid();
        void distributeToLocal();

    private:
        int structuredGridNodeSet(
            const int nxNodes, const int nyNodes, const int nzNodes,
            const int i, const int j, const int k, const int p);
        double structuredGridCoordinateSet(
            const double dx, const double dy, const double dz,
            const int i, const int j, const int k, const int d);
        
};

struct VoxelInfo
{
    std::vector<double> v;
    std::vector<double> center;
    std::vector<int> cellChildren;

    void setNearCell(Grid &grid, double length, int &dim);
    void averageVelocity(Node &node, const int nNodesInCell, const int dim);
    void gaussIntegral(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent, 
                       std::vector<std::vector<double>> &velCurrent, const int &nNodesInCell, 
                       const double &detJ, const double &weight, const double &dim);
};

class ObservedGrid
{
    public:
        ObservedGrid(){}  
        ObservedGrid(Config &conf);

        int nx, ny, nz;
        double dx, dy, dz; 
        double lx, ly, lz;

        int nCellsGlobal;
        int nNodesGlobal;
        int nNodesInCell;
        
        inline VoxelInfo& operator()(int x)
        { return data[x]; }

        inline VoxelInfo& operator()(int y, int x)
        { return data[y * dx + x]; }

        inline VoxelInfo& operator()(int z, int y, int x)
        { return data[z * dx * dy + y * dx + x]; }

        inline int size()
        { return data.size(); }

        inline void resize(int n)
        { data.resize(n); }

        std::vector<VoxelInfo> data;
};

#endif