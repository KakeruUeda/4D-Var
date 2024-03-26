#ifndef GRID_H
#define GRID_H

#include <iostream>
#include "Config.h"
#include "Cell.h"

class Grid
{
    public:
        Grid(Config &conf);
        virtual ~Grid(){};

        GridType gridType;

        Array1D<size_t> nx;
        Array1D<double> lx;
        Array1D<double> dx;

        size_t dim;
        size_t nNodesGlobal;
        size_t nCellsGlobal;

    private:
        Cell<CellInfo> cell;

};

#endif