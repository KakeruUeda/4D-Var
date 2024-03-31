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
        Cell(){};
        Cell(Config &conf) :
        nCellsGlobal(conf.nCellsGlobal), data(conf.nCellsGlobal)
        {  
            for(size_t i=0; i<nCellsGlobal; i++) 
                data[i].nNodesInCell = conf.nNodesInCell;
             for(size_t i=0; i<nCellsGlobal; i++) 
                data[i].node.resize(conf.nNodesInCell);
        }
        ~Cell(){}
        
        inline T& operator()(size_t n)
        { return data[n]; }

        inline void resize(size_t n)
        { data.resize(n); }

        inline void setZero()
        { for(T& n : data) n = 0; }

        size_t nCellsGlobal;
        size_t meshType;

    private:
	    std::vector<T> data;
};


struct CellInfo
{
    public:
        size_t nNodesInCell;
        Array1D<size_t> node;

        inline void setArrayZero(size_t n);
};



#endif