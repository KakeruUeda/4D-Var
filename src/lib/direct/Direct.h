#ifndef DIRECT_H
#define DIRECT_H

#include <iostream>
#include "Config.h"
#include "Cell.h"

struct EstimatedVariable
{
    public:
        Array1D<double> u, v ,w;
};

class Direct
{
    public:
        Direct(Config &conf) : cell(conf){}
        ~Direct(){}

        Cell<CellProperty> cell;
        Array1D<double> u, v ,w;
        double Re, pho, mu, nu;

};

#endif

