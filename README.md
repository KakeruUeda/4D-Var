# 4D-Var
## Overview
This C++ code supports both foward and inverse analyses in CFD using the Message Passing Interface (MPI). The inverse routine addresses forward and adjoint equations to optimize control variables, including the inlet Dirichlet boundary condition and the initial velocity field. These physical equations are discretized on an orthogonal grid employing the Finite Element Method (FEM).
## Dependencies
・PETSc: The Portable, Extensible Toolkit for Scientific Computation <br>
・METIS: Parallel Domain Partitioning <br>
・TextParser: Text Parsing Library for input <br>

See: <br>
https://petsc.org/release/ <br>
https://github.com/KarypisLab/METIS <br>
https://github.com/avr-aics-riken/TextParser <br>

## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> <tp_name>.tp petsc_options.dat
## Applications
・**UnsteadyNavierStokesSolver**: Solves Unsteady Navier-Stokes equation. <br>
・**4DVar (Four-Dimensional Variational Data Assimilation)**: Solves inverse problem. <br> 
## Examples
&nbsp;&nbsp;&nbsp; Foward Solver (8 MPI Processes) <br>
<img src="images/vessel_group.png" alt="Image description" width="600"> <br>
&nbsp;&nbsp;&nbsp; Example - MPI Performance <br>
<img src="images/mpi_performance.png" alt="Image description" width="270"> <br>