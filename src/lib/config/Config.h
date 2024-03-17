#ifndef CONFIGURE_H
#define CONFIGURE_H

#include <iostream>
#include <vector>
#include <cassert>
 
class Config
{
    public:
        Config(){};
        ~Config(){};
        
        size_t nx, ny, nz;
        double lx, ly, lz;
        double dx, dy, dz;

        size_t nCellsGlobal;
        size_t nNodesGlobal;
        size_t nNodesInCellTmp;

        std::vector<double> nNodesInCell;

        void readStructuredParam();
        void readUnstructuredParam();
        void set();
};

#endif