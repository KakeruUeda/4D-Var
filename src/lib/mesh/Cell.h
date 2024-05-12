#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <cassert>
#include "Array.h"
#include "Config.h"
#include "VTKCellType.h"

struct CellInfo
{
    public:
        VTKCellType cellType;
        int nNodesInCell;
        int subId;
        Array1D<int> node;
        Array2D<double> x;

        inline void setArrayZero(int n);
};

class Cell
{
    public:
        Cell(){}
        Cell(Config conf) :
        nCellsGlobal(conf.nCellsGlobal), data(conf.nCellsGlobal){}
        ~Cell(){}
        
        inline CellInfo& operator()(int n)
        { return data[n]; }

        inline void resize(int n)
        { data.resize(n); }

        int nCellsGlobal;

        void initialize(Config conf);

    private:
	    std::vector<CellInfo> data;
};



#endif