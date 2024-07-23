/**
 * @file Postprocess.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Postprocess.h"

void Postprocess::extractOutletVelocity(DirectProblem &direct, std::vector<int> sortNode)
{
    if(mpi.myId > 0) return;

    int nodeCount = 0;
    double xMax = 4e0;
    double eps = 1e-8;
    xMax = xMax - eps;

    std::ofstream outDirichletOutlet(direct.outputDir + "/dat/velocityDirichletPoiseuille.dat");
    for(int in=0; in<direct.grid.node.nNodesGlobal; in++){
        if(direct.grid.node.x[in][0] > xMax){
            int index = sortNode[in] - 128;
            outDirichletOutlet << sortNode[in] << " ";
            for(int d=0; d<direct.dim; d++){
                outDirichletOutlet << direct.grid.node.v[in][d] << " ";
            }
            outDirichletOutlet << std::endl;
        }
    }
    outDirichletOutlet.close();

}

void Postprocess::createData(DirectProblem &direct)
{
    //voxel.range = 5e-1 * voxel.dx;
    voxel.range = 5e-1 * sqrt(voxel.dx * voxel.dx + voxel.dy * voxel.dy + voxel.dz * voxel.dz);
    //voxel.range = 2 * voxel.dx + voxel.dx;
    
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
                    voxel(k, j, i).average(direct.grid.cell, direct.snap.v[t], t, direct.dim);
                }
            }
        }
    }
    
    for(int t=0; t<direct.snap.nSnapShot; t++){
        if(mpi.myId == 0){
            std::string vtiFile;
            vtiFile = direct.outputDir + "/data/data" + to_string(t) + ".vti";
            direct.grid.vtk.exportDataVTI(vtiFile, voxel, t, direct.dim);
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
            std::vector<std::vector<double>> velRef;
            std::vector<std::vector<double>> velRef2;
            int n = (direct.grid.nx+1) * (direct.grid.ny+1) * (direct.grid.nz+1);
            VecTool::resize(velRef, n, direct.dim);
            for(int in=0; in<direct.grid.node.nNodesGlobal; in++){
                velRef[direct.grid.node.sortNode[in]] = direct.snap.v[t][in];
            }
            std::ofstream outReference(direct.outputDir + "/data/reference" + to_string(t) + ".dat");
            for(int k=0; k<direct.grid.nz+1; k++){
                for(int j=0; j<direct.grid.ny+1; j++){
                    for(int i=0; i<direct.grid.nx+1; i++){
                        int index = i + j * (direct.grid.nx+1) + k * (direct.grid.nx+1) * (direct.grid.ny+1);
                        if(i >= direct.grid.nx/2){
                            for(int d=0; d<direct.dim; d++){
                                outReference << velRef[index][d] << " ";
                            }
                            velRef2.push_back(velRef[index]);
                            outReference << std::endl;
                        }
                    }
                }
            }
            outReference.close();
            
            std::string vtiFile;
            vtiFile = direct.outputDir + "/data/reference" + to_string(t) + ".vti";
            direct.grid.vtk.exportNodeVTI(vtiFile, velRef2, direct.grid.nx/2, direct.grid.ny, direct.grid.nz, direct.grid.dx, direct.grid.dy, direct.grid.dz);
        }
    }

}