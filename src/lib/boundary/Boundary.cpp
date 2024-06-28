#include "Boundary.h"

void StructuredBoundaryFace::setNodesOnBoundaryFace(int nxNodes, int nyNodes, int nzNodes)
{
    if(bdFaceStr == "top"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(j == nyNodes-1){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
    if(bdFaceStr == "bottom"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(j == 0){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
    if(bdFaceStr == "left"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(i == 0){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
    if(bdFaceStr == "right"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(i == nxNodes-1){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
    if(bdFaceStr == "front"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(k == 0){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
    if(bdFaceStr == "back"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(k == nzNodes-1){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
}


void StructuredBoundaryFace::setDirichletInfo(std::vector<std::string> bdType, 
                                              std::vector<std::vector<double>> bdValue, 
                                              int dim, int bdIndex)
{
    dirichletType.resize(node.size());
    dirichletValue.resize(node.size());

    for(int i=0; i<dirichletType.size(); i++)
        dirichletType[i] = bdType[bdIndex];

    for(int i=0; i<dirichletValue.size(); i++)
        for(int d=0; d<bdValue[bdIndex].size(); d++)
            dirichletValue[i].push_back(bdValue[bdIndex][d]);
}

void DirichletBoundary::initialize(Config &conf)
{
    vDirichlet = conf.vDirichlet;
    pDirichlet = conf.pDirichlet;
}

void DirichletBoundary::initializeAdjoint(Config &conf)
{
    vDirichlet = conf.vDirichlet;
    pDirichlet = conf.pDirichlet;
    controlBoundaryMap = conf.controlBoundaryMap;
    controlCellMap = conf.controlCellMap;
    controlNodeInCell = conf.controlNodeInCell;
}

void DirichletBoundary::assignDirichletBCs(Node &node, int &dim)
{
    int dofCurrentTmp, dofCurrent;
    int count;

    for(auto &pair : vDirichlet){
        dofCurrentTmp = 0;
        dofCurrent = 0; 
        count = 0;
        for(int i=0; i<pair.first; i++)
            dofCurrentTmp += node.nDofsOnNode[i];
        for(auto &value : pair.second){
            dofCurrent = dofCurrentTmp + count;
            dirichletBCsValue[dofCurrent] = value;
            dirichletBCsValueInit[dofCurrent] = value;
            count++;
        }
    }

    for(auto &pair : vDirichletNew){
        dofCurrentTmp = 0;
        dofCurrent = 0; 
        count = 0;
        for(int i=0; i<pair.first; i++)
            dofCurrentTmp += node.nDofsOnNodeNew[i];
        for(auto &value : pair.second){
            dofCurrent = dofCurrentTmp + count;
            dirichletBCsValueNew[dofCurrent] = value;
            dirichletBCsValueNewInit[dofCurrent] = value;
            count++;
        }
    }

    for(auto &pair : pDirichlet){
        dofCurrentTmp = 0;
        dofCurrent = 0; 
        count = 0;
        for(int i=0; i<pair.first; i++)
            dofCurrentTmp += node.nDofsOnNode[i];
        dofCurrent = dofCurrentTmp + dim;
        dirichletBCsValue[dofCurrent] = pair.second;
        dirichletBCsValueInit[dofCurrent] = pair.second;
    }

    for(auto &pair : pDirichletNew){
        dofCurrentTmp = 0;
        dofCurrent = 0; 
        count = 0;
        for(int i=0; i<pair.first; i++)
            dofCurrentTmp += node.nDofsOnNodeNew[i];
        dofCurrent = dofCurrentTmp + dim;
        dirichletBCsValueNew[dofCurrent] = pair.second;
        dirichletBCsValueNewInit[dofCurrent] = pair.second;
    }
}

void DirichletBoundary::assignPulsatileBCs(const double &t, const double &dt, 
                                           const double &T, const int &nDofsGlobal)
{
    double timeNow = t * dt;
    double pulse = 0.25 * cos((2e0 * PI / T) * timeNow) + 0.75;
    for(int id=0; id<nDofsGlobal; id++){
        if(dirichletBCsValueNew[id] > 0)
            dirichletBCsValueNew[id] = dirichletBCsValueNewInit[id] * pulse;
    }
}

void DirichletBoundary::applyDirichletBCs(Cell &cell, PetscSolver &petsc)
{
    std::vector<int> vecTmp;

    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        if(cell(ic).subId == mpi.myId){
            int nDofsInCell = cell(ic).dofsMap.size();
            PetscScalar  FlocalTmp[nDofsInCell];
            PetscScalar  KlocalTmp[nDofsInCell * nDofsInCell];

            for(int i=0; i<nDofsInCell; i++)  FlocalTmp[i] = 0e0;
            for(int i=0; i<nDofsInCell*nDofsInCell; i++)  KlocalTmp[i] = 0e0;

            vecTmp = cell(ic).dofsMap;
            for(int i=0; i<nDofsInCell; i++){
                if(cell(ic).dofsBCsMap[i] == -1){
                    KlocalTmp[i + i * nDofsInCell] = 1;
                    FlocalTmp[i] = dirichletBCsValueNew[cell(ic).dofsMap[i]];
                }
            }
            MatSetValues(petsc.mtx, nDofsInCell, &vecTmp[0], nDofsInCell, &vecTmp[0], KlocalTmp, INSERT_VALUES);
            VecSetValues(petsc.rhsVec, nDofsInCell, &vecTmp[0], FlocalTmp, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(petsc.mtx, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(petsc.mtx, MAT_FLUSH_ASSEMBLY);

    VecAssemblyBegin(petsc.rhsVec);
    VecAssemblyEnd(petsc.rhsVec);
}
