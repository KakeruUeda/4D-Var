#include "Cell.h"

inline void CellInfo::setArrayZero(int n)
{   
    assert(n != 0);
}

void Cell::initialize(Config conf)
{
    for(int ic=0; ic<nCellsGlobal; ic++) 
        data[ic].nNodesInCell = conf.nNodesInCell;
    for(int ic=0; ic<nCellsGlobal; ic++) 
        data[ic].node.resize(conf.nNodesInCell);

    for(int ic=0; ic<nCellsGlobal; ic++)
        for(int p=0; p<data[ic].node.size(); ic++)
            data[ic].node(p) = conf.cell[ic][p];

    if(conf.nNodesInCell == 4)
        for(int ic=0; ic<nCellsGlobal; ic++)
            data[ic].cellType = VTK_QUAD;
    if(conf.nNodesInCell == 8)
        for(int ic=0; ic<nCellsGlobal; ic++)
            data[ic].cellType = VTK_VOXEL;

}
