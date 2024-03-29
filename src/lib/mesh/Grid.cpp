#include "Grid.h"

Grid::Grid(Config &conf) : 
cell(conf), nx(conf.dim), lx(conf.dim), dx(conf.dim), dim(conf.dim)
{   
    if(conf.gridTypeString == "Structured")
        gridType = GridType::STRUCTURED;
    else if(conf.gridTypeString == "Unstructured")
        gridType = GridType::UNSTRUCTURED;
}

void Grid::initialize(Config &conf)
{    
    switch(gridType)
    {
        case GridType::STRUCTURED:
            initializeStructuredGrid(conf);
            break;
        case GridType::UNSTRUCTURED:
            initializeUnstructuredGrid(conf);
            break;
        default:
            initializeStructuredGrid(conf);
            break;
    }
    
}

void Grid::initializeStructuredGrid(Config &conf)
{
    for(size_t i=0; i<dim; i++)
    {
        nx(i) = conf.nx(i);
        lx(i) = conf.lx(i);
        nCellsGlobal *= nx(i);
        nNodesGlobal *= (nx(i) + 1);
        dx(i) = lx(i) / double(nx(i));
    }
}

void Grid::initializeUnstructuredGrid(Config &conf)
{
    
}

/*
void Grid::initializeWholeDomain(Config &conf)
{

    if(gridType == GridType::STRUCTURED)
    {
        for(size_t i=0; i<dim; i++)
        {
            nx(i) = conf.nx(i);
            lx(i) = conf.lx(i);
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
*/


