#include <iostream>
#include "petscksp.h"
#include "petscmat.h"

enum 
{ 
    SOLVER_EMPTY, 
    PATTERN_OK, 
    INIT_OK, 
    ASSEMBLY_OK, 
    FACTORISE_OK
};

class PetscSolver
{
    public:
        PetscSolver(){};
        virtual ~PetscSolver(){};
        
        PetscErrorCode errpetsc;

        Vec  rhsVec, solnVec, solnVecPrev;
        Mat  mtx; // linear system matrix
        KSP  ksp; // linear solver context
        PC   pc;  // preconditioner context
        
        PetscInt nRow, nCol, nnz;

        int currentStatus;

        int nnz_max_row;
        PetscInt  *diag_nnz, *offdiag_nnz;

        int initialize(int sizeLocal, int sizeGlobal);
};