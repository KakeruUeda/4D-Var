#include "PetscSolver.h"

PetscSolver::~PetscSolver()
{
    //KSPDestroy(&ksp);
    //VecDestroy(&solnVec);
    //VecDestroy(&rhsVec);
    //MatDestroy(&mtx);
}


int PetscSolver::initialize(int sizeLocal, int sizeGlobal)
{
    nRow = nCol = sizeGlobal;
    int dummy = 50;

    // Create PETSc vector
    errpetsc = VecCreate(PETSC_COMM_WORLD, &solnVec);
    CHKERRQ(errpetsc);
    errpetsc = VecSetSizes(solnVec, sizeLocal, sizeGlobal);
    CHKERRQ(errpetsc);
    errpetsc = VecSetFromOptions(solnVec);
    CHKERRQ(errpetsc);
    errpetsc = VecDuplicate(solnVec, &rhsVec);
    CHKERRQ(errpetsc);
    errpetsc = VecSetOption(rhsVec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

    // Create PETSc matrix
    errpetsc = MatCreate(PETSC_COMM_WORLD, &mtx);
    CHKERRQ(errpetsc);
    errpetsc = MatSetSizes(mtx, sizeLocal, sizeLocal, sizeGlobal, sizeGlobal);
    CHKERRQ(errpetsc);

    errpetsc = MatSetFromOptions(mtx);
    CHKERRQ(errpetsc);
    errpetsc = MatMPIAIJSetPreallocation(mtx, nnz_max_row, NULL, nnz_max_row, NULL);
    CHKERRQ(errpetsc);
    errpetsc = MatSeqAIJSetPreallocation(mtx, nnz_max_row, NULL);
    CHKERRQ(errpetsc);
    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRQ(errpetsc);
    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    CHKERRQ(errpetsc);
    errpetsc = MatSetOption(mtx, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    CHKERRQ(errpetsc);


    // Create the KSP context
    errpetsc = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(errpetsc);

    
    // Set the operators for the KSP context
    errpetsc = KSPSetOperators(ksp, mtx, mtx);
    CHKERRQ(errpetsc);

    //  Set whether to use non-zero initial guess or not
    //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    // Set KSP options from the input file
    // This is convenient as it allows to choose different options
    // from the input files instead of recompiling the code
    errpetsc = KSPSetFromOptions(ksp); 
    CHKERRQ(errpetsc);
    errpetsc = KSPGetPC(ksp, &pc);
    CHKERRQ(errpetsc);
    errpetsc = PCSetFromOptions(pc);
    CHKERRQ(errpetsc);
    

    currentStatus = SOLVER_EMPTY;

    return 0;
}

void PetscSolver::setMatAndVecZero(Cell &cell)
{
    int size = 0;
    std::vector<int> vecTmp;

    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        int nDofsInCell = cell(ic).dofsMap.size();
        PetscScalar  FlocalTmp[nDofsInCell];
        PetscScalar  KlocalTmp[nDofsInCell * nDofsInCell];
        for(int i=0; i<nDofsInCell; i++)  FlocalTmp[i] = 0e0;
        for(int i=0; i<nDofsInCell*nDofsInCell; i++)  KlocalTmp[i] = 0e0;
        
        if(cell(ic).subId == mpi.myId){
            size = cell(ic).dofsMap.size();
            vecTmp = cell(ic).dofsMap;
            MatSetValues(mtx, size, &vecTmp[0], size, &vecTmp[0], KlocalTmp, INSERT_VALUES);
            VecSetValues(rhsVec, size, &vecTmp[0], FlocalTmp, INSERT_VALUES);
        }
    }
}

void PetscSolver::initialAssembly()
{
    MatAssemblyBegin(mtx, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mtx, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhsVec);
    VecAssemblyEnd(rhsVec);
}

void PetscSolver::setValueZero()
{
    MatZeroEntries(mtx);
    VecZeroEntries(rhsVec);
    MatAssemblyBegin(mtx, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mtx, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhsVec);
    VecAssemblyEnd(rhsVec);
}

void PetscSolver::setValue(std::vector<int> &lhsRow, std::vector<int> &lhsColumn,
                           std::vector<int> &rhs, MatrixXd &Klocal, VectorXd &Flocal)
{
    int size1 = lhsRow.size();
    int size2 = lhsColumn.size();
    int size3 = rhs.size();
    MatrixXdRM Klocal2 = Klocal;
    VecSetValues(rhsVec, size3, &rhs[0], &Flocal[0], ADD_VALUES);
    MatSetValues(mtx,    size1, &lhsRow[0], size2, &lhsColumn[0], &Klocal2(0, 0), ADD_VALUES);
}

void PetscSolver::setMatValue(std::vector<int> &lhsRow, std::vector<int> &lhsColumn, MatrixXd &Klocal)
{
    int size1 = lhsRow.size();
    int size2 = lhsColumn.size();
    MatrixXdRM Klocal2 = Klocal;
    MatSetValues(mtx, size1, &lhsRow[0], size2, &lhsColumn[0], &Klocal2(0, 0), ADD_VALUES);
}

void PetscSolver::setVecValue(std::vector<int> &rhs, VectorXd &Flocal)
{
    int size = rhs.size();
    VecSetValues(rhsVec, size, &rhs[0], &Flocal[0], ADD_VALUES);
}

int PetscSolver::solve()
{  
    //if (currentStatus != FACTORISE_OK)
    //{
    //  cerr << " SolverPetsc::solve ... factorise matrix first!" << endl;
    //  return -1;
    //}
    errpetsc = MatAssemblyBegin(mtx, MAT_FINAL_ASSEMBLY);
    CHKERRQ(errpetsc);
    errpetsc = MatAssemblyEnd(mtx, MAT_FINAL_ASSEMBLY);
    CHKERRQ(errpetsc);
    //errpetsc = MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);
    //CHKERRQ(errpetsc);
    
    errpetsc = VecAssemblyBegin(rhsVec);
    CHKERRQ(errpetsc);
    errpetsc = VecAssemblyEnd(rhsVec);
    CHKERRQ(errpetsc);
    //errpetsc = VecView(rhsVec,PETSC_VIEWER_STDOUT_WORLD);
    //CHKERRQ(errpetsc);
    
    errpetsc = VecAssemblyBegin(solnVec);
    CHKERRQ(errpetsc);
    errpetsc = VecAssemblyEnd(solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecZeroEntries(solnVec);
    CHKERRQ(errpetsc);
  
    errpetsc = KSPSolve(ksp, rhsVec, solnVec);
    CHKERRQ(errpetsc);
  
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
 
    PetscInt its;
 
    errpetsc = KSPGetIterationNumber(ksp, &its);
    CHKERRQ(errpetsc);
 
    if(reason < 0){
        PetscPrintf(MPI_COMM_WORLD, "\n Divergence... %d iterations. \n", its);
        std::cout <<  reason << std::endl;
        exit(1);
    }else{
        //PetscPrintf(MPI_COMM_WORLD, "\n Convergence in %d iterations. \n", its);
    }

    return 0;
}

double PetscSolver::vectorNorm(const int num)
{
    double norm = 0e0;
    for(int i=0; i<num; i++){
        norm += solution[i] * solution[i];
    }
    return sqrt(norm);
}
