#include "DirectProblem.h"

DirectProblem::DirectProblem(Config conf) : grid(conf)
{
    grid.cell.initialize(conf);
    grid.node.initialize(conf);
    grid.dirichlet.initialize(conf);
}

void DirectProblem::runSimulation()
{
    petsc.initialize();
    
    if(mpi.nId == 1) 
        prepareSerialMatrix();
    else if(mpi.nId > 1)
        prepareParallelMatrix();
    
    return;
}
