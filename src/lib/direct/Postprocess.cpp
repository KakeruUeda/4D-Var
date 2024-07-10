/**
 * @file Postprocess.cpp
 * @author k.ueda
 * @date Jun, 2024
 */

#include "Postprocess.h"

void Postprocess::extractOutletVelocity(DirectProblem &direct)
{
    if(mpi.myId > 0) return;

    int nodeCount = 0;
    double xMax = 2e0;
    double eps = 1e-8;
    xMax = xMax - eps;
    int index = 0;

    std::ofstream outDirichletOutlet(direct.outputDir + "/dat/velocityDirichletPoiseuille.dat");
    for(int in=0; in<direct.grid.node.nNodesGlobal; in++){
        if(direct.grid.node.x[in][0] > xMax){
            outDirichletOutlet << index << " ";
            for(int d=0; d<direct.dim; d++){
                outDirichletOutlet << direct.grid.node.v[in][d] << " ";
            }
            outDirichletOutlet << std::endl;
            index = index + 33;
        }
    }
    outDirichletOutlet.close();

    index = 0;
    std::ofstream outDirichletOutlet2(direct.outputDir + "/dat/velocityDirichletPoiseuille2.dat");
    for(int in=0; in<direct.grid.node.nNodesGlobal; in++){
        if(direct.grid.node.x[in][0] > xMax){
            outDirichletOutlet2 << index << " ";
            for(int d=0; d<direct.dim; d++){
                outDirichletOutlet2 << direct.grid.node.v[in][d] << " ";
            }
            outDirichletOutlet2 << std::endl;
            index = index + 65;
        }
    }
    outDirichletOutlet2.close();

}

void Postprocess::createData(DirectProblem &direct)
{
    voxel.range = 5e-1 * sqrt(voxel.dx * voxel.dx + voxel.dy * voxel.dy + voxel.dz * voxel.dz);
    for(int k=0; k<voxel.nz; k++){
        for(int j=0; j<voxel.ny; j++){
            for(int i=0; i<voxel.nx; i++){
                for(int d=0; d<direct.dim; d++){
                    if(d == 0) voxel(k, j, i).center[d] = voxel.xOrigin + (5e-1 + i) * voxel.dx;
                    if(d == 1) voxel(k, j, i).center[d] = voxel.yOrigin + (5e-1 + j) * voxel.dy;
                    if(d == 2) voxel(k, j, i).center[d] = voxel.zOrigin + (5e-1 + k) * voxel.dz;
                }
                voxel(k, j, i).setNearCell(direct.grid.node, direct.grid.cell, voxel.range, direct.dim);
                for(int t=0; t<direct.snap.nSnapShot; t++){
                    voxel(k, j, i).average(direct.grid.cell, direct.snap.v[t],
                                                   t, direct.grid.cell.nNodesInCell, direct.dim);
                }
            }
        }
    }
    
    for(int t=0; t<direct.snap.nSnapShot; t++){
        if(mpi.myId == 0){
            std::string vtiFile;
            vtiFile = direct.outputDir + "/data/data" + to_string(t) + ".vti";
            direct.grid.output.exportDataVTI(vtiFile, voxel, t, direct.dim);
        }
    }

    if(mpi.myId == 0){
        for(int t=0; t<direct.snap.nSnapShot; t++){
            std::ofstream outData(direct.outputDir + "/data/data" + to_string(t) + ".dat");
            //outData << voxel.nx << " " << voxel.ny << " " << voxel.nz << std::endl;
            //outData << voxel.lx << " " << voxel.ly << " " << voxel.lz << std::endl;
            for(int k=0; k<voxel.nz; k++){
                for(int j=0; j<voxel.ny; j++){
                    for(int i=0; i<voxel.nx; i++){
                        outData << i << " " << j << " " << j << " ";
                        for(int d=0; d<direct.dim; d++){
                            outData << voxel(k, j, i).vCFD[t][d] << " ";
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
                double pointX = direct.grid.node.x[in][0];
                if(pointX > 1e0 - EPS){
                    for(int d=0; d<direct.dim; d++){
                        outReference << direct.snap.v[t][in][d] << " ";
                    }
                    outReference << std::endl;
                }
            }
            outReference.close();
        }
    }

}