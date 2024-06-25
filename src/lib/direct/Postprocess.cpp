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

void Postprocess::createData(DirectProblem &direct)
{
    double length = 5e-1 * sqrt(voxel.dx * voxel.dx + voxel.dy * voxel.dy + voxel.dz * voxel.dz);
    for(int k=0; k<voxel.nz; k++){
        for(int j=0; j<voxel.ny; j++){
            for(int i=0; i<voxel.nx; i++){
                for(int d=0; d<direct.dim; d++){
                    if(d == 0) voxel(k, j, i).center[d] = voxel.xOrigin + (5e-1 + i) * voxel.dx;
                    if(d == 1) voxel(k, j, i).center[d] = voxel.yOrigin + (5e-1 + j) * voxel.dy;
                    if(d == 2) voxel(k, j, i).center[d] = voxel.zOrigin + (5e-1 + k) * voxel.dz;
                }
                voxel(k, j, i).setNearCell(direct.grid.node, direct.grid.cell, length, direct.dim);
                for(int t=0; t<direct.snap.nSnapShot; t++){
                    voxel(k, j, i).averageVelocity(direct.grid.cell, direct.snap.v[t],
                                                   t, direct.grid.cell.nNodesInCell, direct.dim);
                }
            }
        }
    }
    
    for(int t=0; t<direct.snap.nSnapShot; t++){
        if(mpi.myId == 0){
            std::string vtiFile;
            vtiFile = direct.outputDir + "/data/data" + to_string(t) + ".vti";
            direct.grid.vti.exportDataVTI(vtiFile, voxel, t, direct.dim);
        }
    }

    if(mpi.myId == 0){
        for(int t=0; t<direct.snap.nSnapShot; t++){
            std::ofstream outData(direct.outputDir + "/data/data" + to_string(t) + ".dat");
            outData << voxel.nx << " " << voxel.ny << " " << voxel.nz << std::endl;
            outData << voxel.lx << " " << voxel.ly << " " << voxel.lz << std::endl;
            for(int k=0; k<voxel.nz; k++){
                for(int j=0; j<voxel.ny; j++){
                    for(int i=0; i<voxel.nx; i++){
                        outData << i << " " << j << " " << j << " ";
                        for(int d=0; d<direct.dim; d++){
                            outData << voxel(k, j, i).v[t][d] << " ";
                        }
                        outData << std::endl;
                    }
                }
            }
            outData.close();
        }
    }

    if(mpi.myId == 0){
        for(int t=0; t<direct.snap.nSnapShot; t++){
            std::ofstream outReference(direct.outputDir + "/data/reference" + to_string(t) + ".dat");
            for(int in=0; in<direct.grid.node.nNodesGlobal; in++){
                for(int d=0; d<direct.dim; d++){
                    outReference << direct.snap.v[t][in][d] << " ";
                }
                outReference << std::endl;
            }
            outReference.close();
        }
    }

}