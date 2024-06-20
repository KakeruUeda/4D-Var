#include "PetscSolver.h"

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
    errpetsc = KSPSetFromOptions(ksp);    CHKERRQ(errpetsc);

    errpetsc = KSPGetPC(ksp, &pc);    CHKERRQ(errpetsc);

    errpetsc = PCSetFromOptions(pc);    CHKERRQ(errpetsc);

    currentStatus = SOLVER_EMPTY;

    return 0;
}
