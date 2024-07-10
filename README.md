# 4D-Var
## Overview
Fluid Simulation <br>
## Dependencies
・metis <br>
・petsc <br>
・TextParser
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> <tp_name>.tp petsc_options.dat
## Solvers
・strgrid: Create structured grid <br>
・UnsteadyNavierStokesSolver: Solve Unsteady Navier-Stokes <br>
・4DVar: Four-Dimensional Variational Data Assimilation <br>