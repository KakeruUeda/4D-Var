#include "DataGrid.h"

DataGrid::DataGrid(Config &conf) :
nx(conf.nxData), ny(conf.nyData), nz(conf.nzData),
lx(conf.lxData), ly(conf.lyData), lz(conf.lzData),
dx(conf.dxData), dy(conf.dyData), dz(conf.dzData),
nData(conf.nData), nSnapShot(conf.nSnapShot), 
snapInterval(conf.snapInterval),
nCellsGlobal(conf.nCellsDataGlobal), 
nNodesInCell(conf.nNodesInCellData), 
data(conf.nCellsDataGlobal)
{
    if(conf.app == Application::USNS){
        xOrigin = conf.xOrigin; 
        yOrigin = conf.yOrigin; 
        zOrigin = conf.zOrigin;
        for(int ic=0; ic<conf.nCellsDataGlobal; ic++){
            data[ic].vCFD.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].center.resize(conf.dim, 0e0);
        }
    }else if(conf.app == Application::FDVAR){
        for(int ic=0; ic<conf.nCellsDataGlobal; ic++){
            data[ic].vCFD.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].vMRI.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].ve.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].center.resize(conf.dim, 0e0);
        }
    }
}

void DataGrid::initialize(Config &conf, Node &node, Cell &cell, const int &dim)
{   
    for(int t=0; t<conf.nSnapShot; t++){
        for(int k=0; k<nz; k++){
            for(int j=0; j<ny; j++){
                for(int i=0; i<nx; i++){
                    for(int d=0; d<dim; d++){
                        data[k * nx * ny + j * nx + i].vMRI[t][d] 
                        = conf.velocityData[t][k * nx * ny + j * nx + i][d];
                    }
                }
            }
        }
    }

    range = 5e-1 * sqrt(dx*dx + dy*dy + dz*dz);
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                data[k * nx * ny + j * nx + i].center[0] = (5e-1 + i) * dx;
                data[k * nx * ny + j * nx + i].center[1] = (5e-1 + j) * dy;
                data[k * nx * ny + j * nx + i].center[2] = (5e-1 + k) * dz;
                data[k * nx * ny + j * nx + i].setNearCell(node, cell, range, dim);
            }
        }
    }
    
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                data[k * nx * ny + j * nx + i].setCellOnCenterPoint(node, cell, dim);
            }
        }
    }

}

void VoxelInfo::setNearCell(Node &node, Cell &cell, const double &range, const int &dim)
{
    double distance;
    double diff[dim];
    bool flag;

    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        flag = false;
        for(int p=0; p<cell.nNodesInCell; p++){
            distance = 0e0;
            for(int d=0; d<dim; d++){
                diff[d] = node.x[cell(ic).node[p]][d] - center[d];
                distance += diff[d] * diff[d]; 
            } 
            distance = sqrt(distance);
            if(distance < range) flag = true;
        }
        if(flag) cellChildren.push_back(ic);
    }
}

void VoxelInfo::setCellOnCenterPoint(Node &node, Cell &cell, const int &dim)
{
    double distance;
    std::vector<double> diff(dim);
    std::vector<double> minDiff(dim);
    std::vector<double> point(dim);

    double dx, dy, dz;
    
    double minDistance = 1e12;
    isIncluded = true;

    std::vector<std::vector<double>> xCurrent;
    VecTool::resize(xCurrent, cell.nNodesInCell, dim);

    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        for(int p=0; p<cell.nNodesInCell; p++){
            for(int d=0; d<dim; d++){
                xCurrent[p][d] = cell(ic).x[p][d];
            }
        }
        dx = fabs(xCurrent[0][0] - xCurrent[1][0]);
        dy = fabs(xCurrent[0][1] - xCurrent[3][1]);
        dz = fabs(xCurrent[0][2] - xCurrent[4][2]);
        std::vector<double> N(cell.nNodesInCell, 0e0);
        ShapeFunction3D::C3D8_N(N, 0e0, 0e0, 0e0);
        for(int d=0; d<dim; d++){
            point[d] = 0e0;
            for(int p=0; p<cell.nNodesInCell; p++){
                point[d] += N[p] * node.x[cell(ic).node[p]][d];
            }
        }
        distance = 0e0;
        for(int d=0; d<dim; d++){
            diff[d] = fabs(point[d] - center[d]);
            distance += diff[d] * diff[d];
        }
        distance = sqrt(distance);
        if(distance < minDistance){
            minDistance = distance;
            centerCell = ic;
            for(int d=0; d<dim; d++){
                minDiff[d] = diff[d];
            }
        } 
    }
    if(minDiff[0] > dx/2e0) isIncluded = false;
    if(minDiff[1] > dy/2e0) isIncluded = false;
    if(minDiff[2] > dz/2e0) isIncluded = false;

}

void VoxelInfo::interpolate(Node &node, Cell &cell, std::vector<std::vector<double>> &_v, 
                            const int &t, const int &dim)
{
    if(isIncluded == false) return;

    std::vector<std::vector<double>> xCurrent;
    VecTool::resize(xCurrent, cell.nNodesInCell, dim);

    for(int p=0; p<cell.nNodesInCell; p++){
        for(int d=0; d<dim; d++){
            xCurrent[p][d] = cell(centerCell).x[p][d];
        }
    }
    
    double dx = fabs(xCurrent[0][0] - xCurrent[1][0]);
    double dy = dx;
    double dz = dx;
    double point[3];

    std::vector<double> N(cell.nNodesInCell, 0e0);
    ShapeFunction3D::C3D8_N(N, 0e0, 0e0, 0e0);
    for(int d=0; d<dim; d++){
        point[d] = 0e0;
        for(int p=0; p<cell.nNodesInCell; p++){
            point[d] += N[p] * xCurrent[p][d];
        }
    }
        
    double ss = (center[0] - point[0]);
    double tt = (center[1] - point[1]);
    double uu = (center[2] - point[2]);

    ss = ss / (dx / 2e0);
    tt = tt / (dy / 2e0);
    uu = uu / (dz / 2e0);

    if(ss<-1-EPS || ss>1+EPS){
        PetscPrintf(MPI_COMM_WORLD, "\ns interpolation error found.\n");
    }else if(tt<-1-EPS || tt>1+EPS){
        PetscPrintf(MPI_COMM_WORLD, "\nt interpolation error found.\n");
    }else if(uu<-1-EPS || uu>1+EPS){
        PetscPrintf(MPI_COMM_WORLD, "\nu interpolation error found.\n");
    }

    for(int p=0; p<cell.nNodesInCell; p++){
        N[p] = 0e0;
    }
    ShapeFunction3D::C3D8_N(N, ss, tt, uu);

    for(int d=0; d<dim; d++){
        for(int p=0; p<cell.nNodesInCell; p++){
            vCFD[t][d] += N[p] * _v[cell(centerCell).node[p]][d];
        }
    }

}

void VoxelInfo::averageVelocity(Cell &cell, std::vector<std::vector<double>> &_v, 
                                const int &t, const int &nNodesInCell, const int &dim)
{
    if(cellChildren.size() == 0) return;
    
    double weightIntegral = 0e0;
    for(int ic=0; ic<cellChildren.size(); ic++){
        std::vector<std::vector<double>> velCurrent;
        std::vector<std::vector<double>> xCurrent;

        velCurrent.resize(nNodesInCell, std::vector<double>(dim, 0e0));
        xCurrent.resize(nNodesInCell, std::vector<double>(dim, 0e0));

        for(int p=0; p<nNodesInCell; p++){
            for(int d=0; d<dim; d++){
                velCurrent[p][d] = _v[cell(cellChildren[ic]).node[p]][d];
                xCurrent[p][d] = cell(cellChildren[ic]).x[p][d];
            }
        }
        std::vector<double> N;
        std::vector<std::vector<double>> dNdr;
        N.resize(nNodesInCell, 0e0);
        dNdr.resize(nNodesInCell, std::vector<double>(dim, 0e0));
        int nGaussPoint = 2;
        Gauss gauss(nGaussPoint);
        double detJ, weight;
        
        for(int i1=0; i1<nGaussPoint; i1++){
            for(int i2=0; i2<nGaussPoint; i2++){
                for(int i3=0; i3<nGaussPoint; i3++){
                    double dxdr[3][3];
                    ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    MathFEM::calc_dxdr(dxdr, dNdr, xCurrent, nNodesInCell);
                    detJ = MathCommon::calcDeterminant_3x3(dxdr);
                    weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                    gaussIntegral(N, xCurrent, velCurrent, weightIntegral, nNodesInCell, detJ, weight, t, dim);
                }
            }
        }
    }
    for(int d=0; d<dim; d++)
        vCFD[t][d] /= weightIntegral;
}

void VoxelInfo::gaussIntegral(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent, 
                              std::vector<std::vector<double>> &velCurrent, double &weightIntegral, 
                              const int &nNodesInCell, const double &detJ, const double &weight, 
                              const int &t, const int &dim)
{
    for(int d=0; d<dim; d++){
        for(int p=0; p<nNodesInCell; p++){
            vCFD[t][d] += N[p] * velCurrent[p][d] * detJ * weight;
        }
    }
    for(int p=0; p<nNodesInCell; p++)
        weightIntegral += N[p] * detJ * weight;
}
