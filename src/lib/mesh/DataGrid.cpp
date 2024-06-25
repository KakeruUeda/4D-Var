#include "DataGrid.h"

DataGrid::DataGrid(Config &conf) :
nx(conf.nxData), ny(conf.nyData), nz(conf.nzData),
lx(conf.lxData), ly(conf.lyData), lz(conf.lzData),
dx(conf.dxData), dy(conf.dyData), dz(conf.dzData),
nSnapShot(conf.nSnapShot), snapInterval(conf.snapInterval),
nCellsGlobal(conf.nCellsDataGlobal), 
nNodesInCell(conf.nNodesInCellData), 
data(conf.nCellsDataGlobal)
{
    if(conf.app == Application::USNS){
        xOrigin = conf.xOrigin; 
        yOrigin = conf.yOrigin; 
        zOrigin = conf.zOrigin;
        for(int ic=0; ic<conf.nCellsDataGlobal; ic++){
            data[ic].v.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].center.resize(conf.dim, 0e0);
        }
    }else if(conf.app == Application::FDVAR){
        for(int ic=0; ic<conf.nCellsDataGlobal; ic++){
            data[ic].vCFD.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].vMRI.resize(conf.nSnapShot, std::vector<double>(conf.dim, 0e0));
            data[ic].center.resize(conf.dim, 0e0);
        }
    }
}


void VoxelInfo::setNearCell(Node &node, Cell &cell, const double &length, const int &dim)
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
            if(distance < length) flag = true;
        }
        if(flag) cellChildren.push_back(ic);
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
        v[t][d] /= weightIntegral;
}

void VoxelInfo::gaussIntegral(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent, 
                              std::vector<std::vector<double>> &velCurrent, double &weightIntegral, 
                              const int &nNodesInCell, const double &detJ, const double &weight, 
                              const int &t, const int &dim)
{
    for(int d=0; d<dim; d++){
        for(int p=0; p<nNodesInCell; p++){
            v[t][d] += N[p] * velCurrent[p][d] * detJ * weight;
        }
    }
    for(int p=0; p<nNodesInCell; p++)
        weightIntegral += N[p] * detJ * weight;
}
