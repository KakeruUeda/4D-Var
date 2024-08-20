# 4D-Var
## Overview
This C++ code supports both foward and inverse analyses in CFD using the Message Passing Interface (MPI). The inverse routine solves forward and adjoint equations to optimize control variables, including the inlet Dirichlet boundary condition and the initial velocity field. These physical equations are discretized on an orthogonal grid employing the Finite Element Method (FEM).
## Dependencies
・PETSc: Parallel FEM Solver <br>
・METIS: Domain Partitioning for Parallel Computing <br>
・TextParser: Text Parsing Library to Handle Input Parameters <br>

See: <br>
https://petsc.org/release/ <br>
https://github.com/KarypisLab/METIS <br>
https://github.com/avr-aics-riken/TextParser <br>

## Example
    * sh build.sh
    * cd example/inverse/4dvar
    * mpirun -n <num_processes> ./../../../bin/4DVar test.tp petsc_options.dat

## Applications
・**UnsteadyNavierStokesSolver**: Solves Unsteady Navier-Stokes equation. <br>
・**4DVar (Four-Dimensional Variational Data Assimilation)**: Solves inverse problem. <br> 