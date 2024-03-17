#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <cassert>
#include "Array.h"
#include "Config.h"

template <class T>
class Cell
{
    public:
        Cell(Config &conf);
        Cell(){};
        ~Cell(){}
        
        inline T& operator()(size_t n)
        { return data.at(n); }

        size_t nCellsGlobal;
        size_t meshType;

    private:
	    std::vector<T> data;
};


class CellProperty
{
    public:
        size_t nNodesInCell;
        Array1D<int> node;
        Array1D<int> nDofs;
        Array1D<double> x, y, z;
        Array1D<double> u, v, w;

        inline void resizeArrayZero(size_t n);
};


template <class T>
Cell<T>::Cell(Config &conf) : 
    nCellsGlobal(conf.nCellsGlobal), data(conf.nCellsGlobal)
{
    for(size_t i=0; i<nCellsGlobal; i++)
    {
        data.at(i).nNodesInCell = conf.nNodesInCell.at(i);
    }

}

#endif