#include "Grid.h"

Grid::Grid(Config &conf) : 
cell(conf), node(conf), nx(conf.dim), lx(conf.dim), dx(conf.dim), dim(conf.dim)
{   
    if(conf.gridTypeString == "Structured")
        gridType = GridType::STRUCTURED;
    else if(conf.gridTypeString == "Unstructured")
        gridType = GridType::UNSTRUCTURED;
}

void Grid::setStructuredGrid(const size_t nxCells, const size_t nyCells, const size_t nzCells, 
                             const size_t nxNodes, const size_t nyNodes, const size_t nzNodes, 
                             const double dx, const double dy, const double dz,
                             const size_t nNodesInCell, const size_t dim, 
                             Cell<CellInfo> &cell, Node<NodeInfo> &node)
{
    for(size_t k=0; k<nzCells; k++)
        for(size_t j=0; j<nyCells; j++)
            for(size_t i=0; i<nxCells; i++)
                for(size_t p=0; p<nNodesInCell; p++)
                    cell(k * nxCells * nyCells + j * nxCells + i).node(p) 
                        = structuredGridNodeSet(nxNodes, nyNodes, nzNodes, i, j, k, p);

    for(size_t k=0; k<nzNodes; k++)
        for(size_t j=0; j<nyNodes; j++)
            for(size_t i=0; i<nxNodes; i++)
                for(size_t d=0; d<dim; d++)
                    node(k * nxNodes * nyNodes + j * nxNodes + i).x(d) 
                        = structuredGridCoordinateSet(dx, dy, dz, i, j, k, d);

    return;
}

size_t Grid::structuredGridNodeSet(const size_t nxNodes, const size_t nyNodes, const size_t nzNodes,
                                   const size_t i, const size_t j, const size_t k, const size_t p)
{
    std::vector<size_t> nodeSet(8);
    nodeSet.at(0) = i   + j*nxNodes       + k*nxNodes*nyNodes;
    nodeSet.at(1) = i+1 + j*nxNodes       + k*nxNodes*nyNodes;
    nodeSet.at(2) = i+1 + (j+1)*nxNodes + k*nxNodes*nyNodes;
    nodeSet.at(3) = i   + (j+1)*nxNodes + k*nxNodes*nyNodes;
    nodeSet.at(4) = i   + j*nxNodes     + (k+1)*nxNodes*nyNodes;
    nodeSet.at(5) = i+1 + j*nxNodes     + (k+1)*nxNodes*nyNodes;
    nodeSet.at(6) = i+1 + (j+1)*nxNodes + (k+1)*nxNodes*nyNodes;
    nodeSet.at(7) = i   + (j+1)*nxNodes + (k+1)*nxNodes*nyNodes;

    return nodeSet.at(p);
}

double Grid::structuredGridCoordinateSet(const double dx, const double dy, const double dz,
                                         const size_t i, const size_t j, const size_t k, const size_t d)
{
    std::vector<double> coordinateSet(3);
    coordinateSet.at(0) = i * dx;
    coordinateSet.at(1) = j * dy;
    coordinateSet.at(2) = k * dz;

    return coordinateSet.at(d);
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


