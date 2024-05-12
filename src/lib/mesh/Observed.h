#include <iostream>
#include <cassert>
#include "Array.h"
#include "Config.h"

struct ObservedVariable
{
    public:
        Array1D<double> u, v, w;
};

template <class T>
class ObservedGrid
{
    public:
        Observed(Configure &conf) :
        nx(conf.nxObs), ny(conf.nyObs), nz(conf.nzObs), data(conf.nxObs * conf.nyObs * conf.nzObs),
        nCellsGlobal(conf.nzObs * conf.nyObs * conf.nxObs), nNodesInCell(conf.nNodesInCell){}
        Obserbed() :
        nx(0), ny(0), nz(0){}

        inline T& operator()(int n)
        { return data.at(n); }

        inline T& operator()(int i, int, j, int k)
        { return data.at(k * width * height + j * width + i); }

        int nxObs, nyObs, nzObs;
        double dxObs, dyObs, dzObs; 
        double lxObs, lyObs, lzObs;
        int nNodesInCell;

    private:
	    std::vector<T> data;
};
