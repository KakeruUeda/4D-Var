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

enum class GridType
{
    STRUCTURED = 0,
    UNSTRUCTURED = 1
};
 
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

#endif