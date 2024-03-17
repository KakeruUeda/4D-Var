#include "Config.h"

void Config::readStructuredParam()
{
    return;
}

void Config::readUnstructuredParam()
{
    return;
}

void Config::set()
{
    assert(nx != 0 && ny != 0 && nz != 0);
    assert(lx != 0 && ly != 0 && lz != 0);

    dx = lx / (double)nx;
    dy = ly / (double)ny;
    dz = lz / (double)nz;

    nCellsGlobal = nx * ny * nz;
    nNodesGlobal = (nx + 1) * (ny + 1) * (nz + 1);

    nNodesInCell.resize(nCellsGlobal, 0e0);
    for(size_t i=0; i<nCellsGlobal; i++)
        nNodesInCell.at(i) = nNodesInCellTmp;

    return;
}