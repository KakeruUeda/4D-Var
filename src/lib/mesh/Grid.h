#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "Config.h"
#include "Cell.h"
#include "Node.h"
#include "Boundary.h"

enum class GridType
{
    STRUCTURED = 0,
    UNSTRUCTURED = 1
};
 

class Grid
{
    public:
        Grid(){};
        Grid(Config &conf);
        virtual ~Grid(){};

        GridType gridType;
        Cell<CellInfo> cell;
        Node<NodeInfo> node;
        DirichletBoundary boundary;

        size_t nx, ny, nz;
        double lx, ly, lz;
        double dx, dy, dz;

        size_t dim;
        size_t nNodesGlobal;
        size_t nCellsGlobal;

        void setStructuredGrid(
            const size_t nxCells, const size_t nyCells, const size_t nzCells, 
            const size_t nxNodes, const size_t nyNodes, const size_t nzNodes, 
            const double dx, const double dy, const double dz,
            const size_t nNodesInCell, const size_t dim, 
            Cell<CellInfo> &cell, Node<NodeInfo> &node);

    private:
        size_t structuredGridNodeSet(
            const size_t nxNodes, const size_t nyNodes, const size_t nzNodes,
            const size_t i, const size_t j, const size_t k, const size_t p);
        double structuredGridCoordinateSet(
            const double dx, const double dy, const double dz,
            const size_t i, const size_t j, const size_t k, const size_t d);
};

#endif