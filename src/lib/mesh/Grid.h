#ifdef GRID_H
#define GRID_H

#include <iostream>

class StructuredGrid
{
    public:
        StructuredGrid(const ConfigTextParser &params) : 
        nx(params.nx), ny(param.ny), nz(param.nz), 
        lx(params.nx), ly(param.ny), lz(param.nz)
        nNodesInElm(param.nNodesInElm), nCellsGlobal(param.nCellsGlobal)
        { dx = lx / double(nx), dy = ly / double(ny), lz = lz / double(nz) };

        virtual ~StructuredGrid(){};
        
        const size_t nx, ny, nz;
        const double lx, ly, lz;
        const double dx, dy, dz;

        const size_t nNodesGlobal;
        const size_t nCellsGlobal;

};


class UnstructuredGrid
{
};

#endif