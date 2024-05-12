#include <iostream>
#include "petscksp.h"
#include "petscmat.h"

class PetscSolver
{
    public:
        PetscSolver(){};
        virtual ~PetscSolver(){};
        
        PetscErrorCode errpetsc;

        Vec  rhs;
        Mat  mtx; // linear system matrix
        KSP  ksp; // linear solver context
        PC   pc;  // preconditioner context
        
        PetscInt nRow, nCol, nnz;

        int nnz_max_row;
        PetscInt  *diag_nnz, *offdiag_nnz;

        void initialize();
};