#include "Postprocess.h"

void Postprocess::extractOutletVelocity(DirectProblem &direct)
{
    if(mpi.myId > 0) return;

    int nodeCount = 0;
    double xMax = 2e0;
    double eps = 1e-8;
    xMax = xMax - eps;

    std::ofstream outDirichletOutlet(direct.outputDir + "/dat/velocityDirichletPoiseuille.dat");
    for(int in=0; in<direct.grid.node.nNodesGlobal; in++){
        if(direct.grid.node.x[in][0] > xMax){
            int nodeRight = nodeCount - 64;
            outDirichletOutlet << nodeRight << " ";
            for(int d=0; d<direct.dim; d++){
                outDirichletOutlet << direct.grid.node.v[in][d] << " ";
            }
            outDirichletOutlet << std::endl;
        }
        nodeCount++;
    }
}

void Postprocess::makeObservedData(DirectProblem &direct)
{
    double length = 5e-1 * sqrt(obs.dx * obs.dx + obs.dy * obs.dy + obs.dz * obs.dz);
    
    for(int k=0; k<obs.nz; k++){
        for(int j=0; j<obs.ny; j++){
            for(int i=0; i<obs.nx; i++){
                for(int d=0; d<direct.dim; d++){
                    if(d = 0) obs(k, j, i).center[d] = 5e-1 + i * obs.dx;
                    if(d = 1) obs(k, j, i).center[d] = 5e-1 + j * obs.dy;
                    if(d = 2) obs(k, j, i).center[d] = 5e-1 + k * obs.dz;
                }
                obs(k, j, i).setNearCell(direct.grid, length, direct.dim);
                obs(k, j, i).averageVelocity(direct.grid.node, direct.grid.cell.nNodesInCell, direct.dim);
            }
        }
    }

}