# 4D-Var
## Overview
This C++ code supports both direct and inverse analyses in CFD.<br>
## Dependencies
・METIS: Parallel Domain Partitioning <br>
・PETSc: The Portable, Extensible Toolkit for Scientific Computation <br>
・TextParser: Text Parsing Library for input <br>
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> <tp_name>.tp petsc_options.dat
## Applications
・**Strgrid**: Create structured grid <br>
・**UnsteadyNavierStokesSolver**: Solve Unsteady Navier-Stokes <br>
&nbsp;&nbsp;&nbsp; Example - 8 MPI Processes <br>
<img src="images/vessel_group.png" alt="Image description" width="600"> <br>
&nbsp;&nbsp;&nbsp; Example - MPI Performance <br>
<img src="images/mpi_performance.png" alt="Image description" width="270"> <br>
・**4DVar**: Four-Dimensional Variational Data Assimilation <br>
・**FlowRate, MAE**: 4DVar postproces <br>