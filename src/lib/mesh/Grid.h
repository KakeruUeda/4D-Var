#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "Config.h"
#include "Cell.h"
#include "Boundary.h"

enum class GridType
{
    STRUCTURED = 0,
    UNSTRUCTURED = 1
};
 

class Grid
{
    public:
        Grid(Config &conf);
        virtual ~Grid(){};

        GridType gridType;
        Cell<CellInfo> cell;
        Boundary boundary;

        Array1D<size_t> nx;
        Array1D<double> lx;
        Array1D<double> dx;

        size_t dim;
        size_t nNodesGlobal;
        size_t nCellsGlobal;

        void initialize(Config &conf);

    private:
        void initializeStructuredGrid(Config &conf);
        void initializeUnstructuredGrid(Config &conf);

};

#endif