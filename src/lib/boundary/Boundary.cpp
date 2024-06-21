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
    for(int ib=0; ib<conf.vDirichletNode.size(); ib++)
        velocity[ib].node = conf.vDirichletNode[ib];

    for(int ib=0; ib<conf.vDirichletValue.size(); ib++){
        velocity[ib].value.resize(conf.vDirichletValue[ib].size());
        for(int d=0; d<conf.vDirichletValue[ib].size(); d++){
            velocity[ib].value[d] = conf.vDirichletValue[ib][d];
        }
    }

    for(int ib=0; ib<conf.pDirichletNode.size(); ib++)
        pressure[ib].node = conf.pDirichletNode[ib];

    for(int ib=0; ib<conf.pDirichletValue.size(); ib++)
        pressure[ib].value = conf.pDirichletValue[ib];
       
}

void DirichletBoundary::assignDirichletBCs(Node &node, int &dim)
{
    int dofCurrentTmp, dofCurrent;

    // velocity
    for(int ib=0; ib<nNodesVelocity; ib++){
        dofCurrentTmp = 0;
        dofCurrent = 0;
        for(int i=0; i<velocity[ib].node; i++)
            dofCurrentTmp += node.nDofsOnNode[i];
        for(int d=0; d<dim; d++){
            dofCurrent = dofCurrentTmp + d;
            dirichletBCsValue[dofCurrent] = velocity[ib].value[d];
        }
    }

    // velocity new
    for(int ib=0; ib<nNodesVelocity; ib++){
        dofCurrentTmp = 0;
        dofCurrent = 0;
        for(int i=0; i<velocity[ib].nodeNew; i++)
            dofCurrentTmp += node.nDofsOnNode[i];
        for(int d=0; d<dim; d++){
            dofCurrent = dofCurrentTmp + d;
            dirichletBCsValueNew[dofCurrent] = velocity[ib].value[d];
        }
    }

    // prresure
    for(int ib=0; ib<nNodesPressure; ib++){
        dofCurrentTmp = 0;
        dofCurrent = 0;
        for(int i=0; i<pressure[ib].node; i++)
            dofCurrent += node.nDofsOnNode[i];
        for(int d=0; d<dim+1; d++){
            dofCurrent = dofCurrentTmp + d;
            dirichletBCsValue[dofCurrent] = pressure[ib].value;
        }
    }

    // prresure new
    for(int ib=0; ib<nNodesPressure; ib++){
        dofCurrentTmp = 0;
        dofCurrent = 0;
        for(int i=0; i<pressure[ib].nodeNew; i++)
            dofCurrent += node.nDofsOnNode[i];
        for(int d=0; d<dim+1; d++){
            dofCurrent = dofCurrentTmp + d;
            dirichletBCsValueNew[dofCurrent] = pressure[ib].value;
        }
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
