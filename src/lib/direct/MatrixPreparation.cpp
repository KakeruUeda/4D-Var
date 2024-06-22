#include "DirectProblem.h"

void DirectProblem::preprocess()
{
    for(auto &pair : grid.dirichlet.vDirichlet){
        int count = 0;
        for(auto &value : pair.second){
            grid.node.isDirichlet[pair.first][count] = true;
            count++;
        }
    }

    for(auto &pair : grid.dirichlet.pDirichlet)
        grid.node.isDirichlet[pair.first][dim] = true;

    for(int in=0; in<grid.node.nNodesGlobal; in++)
        for(int id=0; id<grid.node.nDofsOnNode[in]; id++)
            if(grid.node.isDirichlet[in][id])
                grid.node.dofsBCsMap[in][id] = -1;

    for(int in=0; in<grid.node.nNodesGlobal; in++)
        for(int id=0; id<grid.node.nDofsOnNode[in]; id++)
                grid.nDofsGlobal++;

    // debug
    if(mpi.myId == 0){
        std::ofstream outIsDirichlet(outputDir + "/dat/isDirichlet.dat");
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int id=0; id<grid.node.nDofsOnNode[in]; id++){
                outIsDirichlet << grid.node.isDirichlet[in][id] << " ";
            }
            outIsDirichlet << std::endl;
        }
        outIsDirichlet.close();
    }

    //debug
    if(mpi.myId == 0){
        std::ofstream outDofsBCsMap(outputDir + "/dat/dofsBCsMap.dat");
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int id=0; id<grid.node.nDofsOnNode[in]; id++){
                outDofsBCsMap << grid.node.dofsBCsMap[in][id] << " ";
            }
            outDofsBCsMap << std::endl;
        }
        outDofsBCsMap.close();
    }

    if(mpi.nId == 1){
        prepareSerialMatrix();
    }else if(mpi.nId > 1){
        grid.divideWholeGrid();
        grid.distributeToLocal();
    }

    // debug
    if(mpi.myId == 0){
        std::ofstream outIsDirichletNew(outputDir + "/dat/isDirichletNew.dat");
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int id=0; id<grid.node.nDofsOnNode[in]; id++){
                outIsDirichletNew << grid.node.isDirichletNew[in][id] << " ";
            }
            outIsDirichletNew << std::endl;
        }
        outIsDirichletNew.close();
    }

    //debug
    if(mpi.myId == 0){
        std::ofstream outDofsBCsMapNew(outputDir + "/dat/dofsBCsMapNew.dat");
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int id=0; id<grid.node.nDofsOnNode[in]; id++){
                outDofsBCsMapNew << grid.node.dofsBCsMapNew[in][id] << " ";
            }
            outDofsBCsMapNew << std::endl;
        }
        outDofsBCsMapNew.close();
    }

    int size = 0;
    for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
        if(grid.cell(ic).subId == mpi.myId){
            size = 0;
            for(int p=0; p<grid.cell.nNodesInCell; p++){
                size += grid.node.nDofsOnNode[grid.cell(ic).nodeNew[p]];
            }
            grid.cell(ic).dofsMap.resize(size);
            grid.cell(ic).dofsBCsMap.resize(size);

            for(int p=0; p<grid.cell.nNodesInCell; p++){
                int i = grid.node.nDofsOnNode[grid.cell(ic).nodeNew[p]] * p;
                int j = grid.cell(ic).nodeNew[p];
                for(int q=0; q<grid.node.nDofsOnNode[grid.cell(ic).nodeNew[p]]; q++){
                    grid.cell(ic).dofsMap[i+q] = grid.node.dofsMapNew[j][q];
                    grid.cell(ic).dofsBCsMap[i+q] = grid.node.dofsBCsMapNew[j][q];
                }
            }
        }
    }

    // debug
    if(mpi.myId == 0){
        std::ofstream outCellDofsMap(outputDir + "/dat/cellDofsBCsMap.dat");
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                for(int p=0; p<grid.cell.nNodesInCell; p++){
                    int i = grid.node.nDofsOnNode[grid.cell(ic).nodeNew[p]] * p;
                    for(int q=0; q<grid.node.nDofsOnNode[grid.cell(ic).nodeNew[p]]; q++){
                        outCellDofsMap << grid.cell(ic).dofsBCsMap[i+q] << " ";
                    }
                    outCellDofsMap << "  ";
                }
                outCellDofsMap << std::endl;
            }
        }
        outCellDofsMap.close();
    }

    int *tt, tmpInt;
    int r, kk, nSize;
    int countDiag, countOffDiag;

    std::vector<std::set<int>> forAssyMatFluid;
    std::set<int>::iterator it;

    forAssyMatFluid.resize(grid.nDofsGlobal);

    for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
        if(grid.cell(ic).subId == mpi.myId){
            tt = &(grid.cell(ic).dofsMap[0]);
            nSize = grid.cell(ic).dofsMap.size();

            for(int i=0; i<nSize; i++){
                r = tt[i];
                if(r != -1){
                    if(r >= grid.rowStart && r <= grid.rowEnd){
                        for(int j=0; j<nSize; j++){
                            if(tt[j] != -1){
                                forAssyMatFluid[r].insert(tt[j]);
                            }
                        }
                    }
                }
            }
        }
    }
  
    PetscMalloc1(grid.nDofsLocal,  &petsc.diag_nnz);
    PetscMalloc1(grid.nDofsLocal,  &petsc.offdiag_nnz);

    kk = 0;
    petsc.nnz_max_row = 0;
    for(int i=grid.rowStart; i<=grid.rowEnd; i++){
        nSize = forAssyMatFluid[i].size();
        petsc.nnz_max_row = std::max(petsc.nnz_max_row, nSize);

        countDiag=0, countOffDiag=0;
        for(it=forAssyMatFluid[i].begin(); it!=forAssyMatFluid[i].end(); it++){   
            tmpInt = *it;
            if(tmpInt >= grid.rowStart && tmpInt <= grid.rowEnd)
                countDiag++;
            else
                countOffDiag++;
        }
        petsc.diag_nnz[kk]    = countDiag;
        petsc.offdiag_nnz[kk] = countOffDiag;
        kk++;
    }

    petsc.initialize(grid.nDofsLocal, grid.nDofsGlobal);

}

void DirectProblem::prepareSerialMatrix()
{

}
