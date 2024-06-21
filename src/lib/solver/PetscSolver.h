#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

#include <iostream>
#include "Eigen.h"
#include "Cell.h"
#include "petscksp.h"
#include "petscmat.h"

enum 
{ 
    SOLVER_EMPTY, PATTERN_OK, 
    INIT_OK, ASSEMBLY_OK, FACTORISE_OK
};

class PetscSolver
{
    public:
        PetscSolver(){};
        virtual ~PetscSolver();
        
        PetscErrorCode errpetsc;

        Vec  rhsVec, solnVec, solnVecPrev;
        Mat  mtx; // linear system matrix
        KSP  ksp; // linear solver context
        PC   pc;  // preconditioner context
        
        PetscInt nRow, nCol, nnz;

        int currentStatus;

        int nnz_max_row;
        PetscInt  *diag_nnz, *offdiag_nnz;

        std::vector<double> solution;

        int initialize(int sizeLocal, int sizeGlobal);
        void setMatAndVecZero(Cell &cell);
        void setValueZero();
        void initialAssembly();
        void setValue(std::vector<int> dofsBCsMap, std::vector<int> dofsMap, 
                      MatrixXd& Klocal, VectorXd& Flocal);
        int solve();
};

#endif