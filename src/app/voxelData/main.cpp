/**
 * @file main.cpp
 * @author K.Ueda
 * @date July, 2024
*/

#include <unistd.h>
#include "DirectProblem.h"
#include "MyMPI.h"
MyMPI mpi;

void createData(Config &conf, Cell &cell, Node &node, DataGrid &voxel, SnapShot &snap);
void createReferenceDA(Config &conf, std::vector<std::vector<double>> &vt, 
                       std::vector<std::vector<double>> &velRef);

int main(int argc, char *argv[])
{
    std::string inputFile = argv[1];
    std::string appName = "VOXELDATA";
    Config conf(inputFile, appName);
    if(conf.isReadingError) return EXIT_FAILURE;

    std::string dir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    conf.outputDir = "output/" + conf.outputDir;
    mkdir(conf.outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = conf.outputDir + "/input";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = conf.outputDir + "/debug";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    Cell cell; Node node;
    SnapShot snap(conf);
    DataGrid voxel(conf);
    
    std::vector<std::vector<std::vector<double>>> vOrig;
    std::vector<std::vector<std::vector<double>>> vRef;
    
    VecTool::resize(vRef, conf.stepMax, conf.nNodesOptGlobal, conf.dim);
    VecTool::resize(snap.v, snap.nSnapShot, conf.nNodesGlobal, conf.dim);
    
    cell.resize(conf.nCellsGlobal);
    node.x.resize(conf.nNodesGlobal, std::vector<double>(conf.dim));

    for(int ic=0; ic<conf.nCellsGlobal; ic++){
        cell(ic).node.resize(conf.nNodesInCell);
        for(int p=0; p<conf.nNodesInCell; p++){
            cell(ic).node[p] = conf.cell[ic][p];
        }
    }

    for(int ic=0; ic<conf.nCellsGlobal; ic++){
        cell(ic).x.resize(conf.nNodesInCell, std::vector<double>(conf.dim));
    }

    for(int ic=0; ic<conf.nCellsGlobal; ic++){
        for(int p=0; p<conf.nNodesInCell; p++){
            for(int d=0; d<conf.dim; d++){
                cell(ic).x[p][d] = conf.node[cell(ic).node[p]][d];
            }
        }
    }
    
    for(int in=0; in<conf.nNodesGlobal; in++){
        for(int d=0; d<conf.dim; d++){
            node.x[in][d] = conf.node[in][d];
        } 
    }

    vOrig.resize(conf.stepMax);
    for(int step=0; step<conf.stepMax; step++){
        std::string velFile = conf.inputDir + "/velocity_" + to_string(step) + ".bin";
        try{
            BIN::importVectorDataBIN(velFile, vOrig[step]);
        }catch(const std::runtime_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }   
    }

    std::cout << vOrig[0].size() << " " << conf.nNodesGlobal << std::endl;

    // debug
    for(int step=0; step<conf.stepMax; step++){
        std::string vtiFile = conf.outputDir + "/debug/velocityOriginal_" + to_string(step) + ".vti";
        VTK::exportVectorPointDataVTI(vtiFile, "velOrig", vOrig[step], conf.nx, conf.ny, conf.nz, conf.dx, conf.dy, conf.dz);
    }

    int snapCount = 0;
    for(int step=0; step<conf.stepMax; step++){
        if(step >= conf.snapTimeBeginItr && (snapCount < conf.nSnapShot)){
            if((step - conf.snapTimeBeginItr) % conf.snapInterval == 0){
                snap.takeSnapShot(vOrig[step], snapCount, conf.nNodesGlobal, conf.dim);
                snapCount++;
            }
        }   
    }

    // debug
    for(int step=0; step<conf.stepMax; step++){
        std::string vtiFile = conf.outputDir + "/debug/velocityOriginal_" + to_string(step) + ".vti";
        VTK::exportVectorPointDataVTI(vtiFile, "velOrig", vOrig[step], conf.nx, conf.ny, conf.nz, conf.dx, conf.dy, conf.dz);
    }

    for(int step=0; step<conf.stepMax; step++){
        createReferenceDA(conf, vOrig[step], vRef[step]);
    }

    // debug
    for(int step=0; step<conf.stepMax; step++){
        std::string vtiFile = conf.outputDir + "/debug/velocityReference_" + to_string(step) + ".vti";
        VTK::exportVectorPointDataVTI(vtiFile, "velRef", vRef[step], conf.nx, conf.ny, conf.nz, conf.dx, conf.dy, conf.dz);
    }

    createData(conf, cell, node, voxel, snap);

    // debug
    for(int step=0; step<snap.nSnapShot; step++){
        std::string vtiFile = conf.outputDir + "/debug/snapShot_" + to_string(step) + ".vti";
        VTK::exportVelocityDataVTI(vtiFile, voxel, step);
    }

    std::cout << "Terminated." << std::endl;

    return EXIT_SUCCESS;
}

void createData(Config &conf, Cell &cell, Node &node, DataGrid &voxel, SnapShot &snap)
{
    voxel.range = 5e-1 * sqrt(voxel.dx * voxel.dx + voxel.dy * voxel.dy + voxel.dz * voxel.dz);
    
    for(int k=0; k<voxel.nz; k++){
        for(int j=0; j<voxel.ny; j++){
            for(int i=0; i<voxel.nx; i++){
                for(int d=0; d<conf.dim; d++){
                    if(d == 0) voxel(k, j, i).center[d] = conf.xOrigin + (5e-1 + i) * voxel.dx;
                    if(d == 1) voxel(k, j, i).center[d] = conf.yOrigin + (5e-1 + j) * voxel.dy;
                    if(d == 2) voxel(k, j, i).center[d] = conf.zOrigin + (5e-1 + k) * voxel.dz;
                }
                voxel(k, j, i).setNearCell(node, cell, voxel.range, conf.dim);
                for(int t=0; t<conf.nSnapShot; t++){
                    voxel(k, j, i).average(cell, snap.v[t], t, conf.dim);
                }
            }
        }
    }
}

void createReferenceDA(Config &conf, std::vector<std::vector<double>> &vt, 
                       std::vector<std::vector<double>> &vRef)
{
    double px, py, pz;
    int ix, iy, iz;
    double s, t, u;
    
    for(int k=0; k<conf.nzOpt+1; k++){
        for(int j=0; j<conf.nyOpt+1; j++){
            for(int i=0; i<conf.nxOpt+1; i++){
                px = conf.xOrigin + i * conf.dxOpt;
                py = conf.yOrigin + j * conf.dyOpt;
                pz = conf.zOrigin + k * conf.dzOpt;
    
                ix = px / conf.nx;
                iy = py / conf.ny;
                iz = pz / conf.nz;
    
                s = ix * conf.dx + 5e-1 * conf.dx;
                t = iy * conf.dy + 5e-1 * conf.dy;
                u = iz * conf.dz + 5e-1 * conf.dz;
                
                s = s / (conf.dx / 2e0);
                t = t / (conf.dy / 2e0);
                u = u / (conf.dz / 2e0);
                
                if(s<-1-EPS || s>1+EPS){
                    std::cout << "s interpolation error." << std::endl;
                }else if(t<-1-EPS || t>1+EPS){
                    std::cout << "t interpolation error." << std::endl;
                }else if(u<-1-EPS || u>1+EPS){
                    std::cout << "u interpolation error." << std::endl;
                }
                
                int n = i + j * (conf.nxOpt+1) + k * (conf.nxOpt+1) * (conf.nyOpt+1);
                int elmcfd = ix + iy * conf.nx + iz * conf.nx * conf.ny;
                
                std::vector<double> N;
                VecTool::resize(N, conf.nNodesInCell);

                ShapeFunction3D::C3D8_N(N, s, t, u);

                for(int d=0; d<conf.dim; d++){
                    for(int p=0; p<conf.nNodesInCell; p++){
                        vRef[n][d] += N[p] * vt[conf.cell[elmcfd][p]][d];
                    }
                }
            }
        }
    }

}
