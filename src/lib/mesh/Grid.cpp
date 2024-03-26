#include "Grid.h"

Grid::Grid(Config &conf) : 
cell(conf), nx(conf.dim), lx(conf.dim), dx(conf.dim),
dim(conf.dim), gridType(conf.gridType)
{
    if(gridType == GridType::STRUCTURED)
    {
        for(size_t i=0; i<dim; i++)
        {
            nx(i) = conf.nx(i);
            lx(i) = conf.lx(i);
        }
        for(size_t i=0; i<dim; i++) 
        {
            nCellsGlobal *= nx(i);
            nNodesGlobal *= (nx(i) + 1);
            dx(i) = lx(i) / double(nx(i));
        }
    }
    else if(gridType == GridType::UNSTRUCTURED)
    {
        nCellsGlobal = conf.nCellsGlobal;
        nNodesGlobal = conf.nNodesGlobal;
    }
    
}
