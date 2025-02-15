/**
 * @file Output.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "InverseProblem.h"

/**
 * @brief Output velocity and pressure.
 */
void InverseProblem::outputFowardSolutions(const int loop)
{
  if(mpi.myId != 0) return;

  if(main.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.updateVelocityVTI(t);
      main.outputVelocityVTI("main", t, loop);
      //main.updatePressureVTI(t);
      //main.outputPressureVTI("main", t, loop);
    }
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.outputVelocityVTU("main", t, loop);
      main.outputPressureVTU("main", t, loop);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }
}

/**
 * @brief Output w, q, l.
 */
void InverseProblem::outputAdjointSolutions(const int loop)
{
  if(mpi.myId != 0) { return; }

  if(adjoint.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < adjoint.timeMax; t++) {
      adjoint.updateSolutionsVTI(t);
      adjoint.outputSolutionsVTI("adjoint", t, loop);
    }
  } else if(adjoint.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < adjoint.timeMax; t++) {
      adjoint.outputSolutionsVTU("adjoint", t, loop);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }
}

/**
 * @brief Output gradients.
 */
void InverseProblem::outputGradients(const int loop)
{
  if(mpi.myId != 0) { return; }

  std::string vtuFile;
  for(int t = 0; t < main.timeMax; t++) {
    vtuFile = main.outputDir + "/other/gradX_" + to_string(loop) + "_" + to_string(t) + ".vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "gradX", main.grid.node, main.grid.cell, gradX, t);
  }
  vtuFile = main.outputDir + "/other/gradX0_" + to_string(loop) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "gradX0", main.grid.node, main.grid.cell, gradX0);
}

/**
 * @brief Output control variables X and X0.
 */
void InverseProblem::outputControlVariables(const int loop)
{
  if(mpi.myId != 0) { return; }

  std::string vtuFile;
  for(int t = 0; t < main.timeMax; t++) {
    vtuFile = main.outputDir + "/other/inletVelocity_" + to_string(loop) + "_" + to_string(t) + ".vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "inletVelocity", main.grid.node, main.grid.cell, X, t);
  }
  vtuFile = main.outputDir + "/other/velocity_initial_" + to_string(loop) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "velocity_initial", main.grid.node, main.grid.cell, X0);
}

/**
 * @brief Update control variables for VTI.
 */
void InverseProblem::updateControlVariablesVTI()
{
  for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      X0vti(main.grid.vecFluidUniqueNodes[in], d) = X0(in, d);
    }
  }
}

/**
 * @brief Output velocity data (vCFD, vMRI, ve).
 */
void InverseProblem::outputVelocityData(const int loop)
{
  if(mpi.myId > 0) return;

  for(int t = 0; t < data.n_mri_step; t++) {
    std::string vtiFile = main.outputDir + "/data/data_" + to_string(loop) + "_" + to_string(t) + ".vti";
    data.exportVTI(vtiFile, t);
  }
}

/**
 * @brief Output time interpolated feedback force.
 */
void InverseProblem::outputFeedbackForce(const int loop)
{
  if(mpi.myId > 0) return;

  std::string vtuFile;
  for(int t = 0; t < main.timeMax; t++) {
    vtuFile = main.outputDir + "/other/feedbackForce" + to_string(loop) + "_" + to_string(t) + ".vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "feedbackForce", main.grid.node, main.grid.cell, adjoint.feedbackForce,
                                     t);
  }
}

/**
 * @brief Output velocity to binary file.
 */
void InverseProblem::outputVelocityBIN(const int loop)
{
  if(mpi.myId > 0) return;

  std::string binFile;
  for(int t = 0; t < main.timeMax; t++) {
    main.updateVelocityVTI(t);
    binFile = main.outputDir + "/bin/velocity_" + to_string(loop) + "_" + to_string(t) + ".bin";
    main.vvti.exportBIN(binFile);
  }

  updateControlVariablesVTI();

  binFile = main.outputDir + "/bin/velocity_initial_" + to_string(loop) + ".bin";
  X0vti.exportBIN(binFile);
}

/**
 * @brief Output optimized variables.
 */
void InverseProblem::outputOptimizedVariables()
{
  if(mpi.myId > 0) return;

  // Optimized velocity field
  if(main.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.updateVelocityVTI(t);
      main.outputVelocityVTI("optimized", t);
      main.outputVelocityVTU("optimized", t);
      // main.updatePressureVTI(t);
      // main.outputPressureVTI("optimized", t);
    }
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.outputVelocityVTU("optimized", t);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }

  // ------ save voxel velocity data ------------------------
  for(int t = 0; t < data.n_mri_step; t++) {
    std::string vtiFile = main.outputDir + "/optimized/data_" + to_string(t) + ".vti";
    data.exportVTI(vtiFile, t);
  }

  for(int t = 0; t < data.n_mri_step; t++) {
    std::string datFile;

    datFile = main.outputDir + "/optimized/v_cfd_" + std::to_string(t) + ".dat";
    data.exportVelCFD_DAT(datFile, t);
    datFile = main.outputDir + "/optimized/v_mri_" + std::to_string(t) + ".dat";
    data.exportVelMRI_DAT(datFile, t);
    datFile = main.outputDir + "/optimized/v_err_" + std::to_string(t) + ".dat";
    data.exportVelError_DAT(datFile, t);
  }
  // ------------------------------------

  // ------ save voxel velocity data (only fluid domain) ------------
  for(int t = 0; t < data.n_mri_step; t++) {
    for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
      if(data.voxel(iv).mask < 0.5) {
        for(int d = 0; d < 3; d++) {
          data.voxel(iv).v_cfd(t, d) = 0e0;
          data.voxel(iv).v_mri(t, d) = 0e0;
          data.voxel(iv).v_err(t, d) = 0e0;
        }
      }
    }
  }

  for(int t = 0; t < data.n_mri_step; t++) {
    std::string vtiFile = main.outputDir + "/optimized/data_fluid_" + to_string(t) + ".vti";
    data.exportVTI(vtiFile, t);
  }

  for(int t = 0; t < data.n_mri_step; t++) {
    std::string datFile;

    datFile = main.outputDir + "/optimized/v_cfd_fluid_" + std::to_string(t) + ".dat";
    data.exportVelCFD_DAT(datFile, t);
    datFile = main.outputDir + "/optimized/v_mri_fluid_" + std::to_string(t) + ".dat";
    data.exportVelMRI_DAT(datFile, t);
    datFile = main.outputDir + "/optimized/v_err_fluid_" + std::to_string(t) + ".dat";
    data.exportVelError_DAT(datFile, t);
  }
  // ------------------------------------

  // To binary file
  if(main.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.updateVelocityVTI(t);
      std::string binFile = main.outputDir + "/optimized/velocity_" + to_string(t) + ".bin";
      main.vvti.exportBIN(binFile);
      binFile = main.outputDir + "/optimized/velocity_fluid_" + to_string(t) + ".bin";
      main.vt.exportBIN(binFile, t);
      
      // ------ TMP ---------
      binFile = main.outputDir + "/optimized/pressure_fluid_" + to_string(t) + ".bin";
      main.pt.exportBIN(binFile, t);
    }
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      std::string binFile = main.outputDir + "/optimized/velocity_" + to_string(t) + ".bin";
      main.vt.exportBIN(binFile, t);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }

  X0vti.allocate(main.grid.node.nNodesStrGlobal, 3);
  
  // Optimized initial velocity field
  if(main.grid.gridType == GridType::STRUCTURED) {
    std::string vtiFile;
    updateControlVariablesVTI();
    vtiFile = main.outputDir + "/optimized/velocity_initial.vti";
    EXPORT::exportVectorPointDataVTI(vtiFile, "velocity_initial", X0vti, main.grid.nx, main.grid.ny, main.grid.nz,
                                     main.grid.dx, main.grid.dy, main.grid.dz);
    std::string vtuFile;
    vtuFile = main.outputDir + "/optimized/velocity_initial.vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "velocity_initial", main.grid.node, main.grid.cell, X0);
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    std::string vtuFile;
    vtuFile = main.outputDir + "/optimized/velocity_initial.vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "velocity_initial", main.grid.node, main.grid.cell, X0);
  } else {
    throw std::runtime_error("Undefined gridType");
  }

  // TO binary file
  if(main.grid.gridType == GridType::STRUCTURED) {
    updateControlVariablesVTI();
    std::string binFile = main.outputDir + "/optimized/velocity_initial.bin";
    X0vti.exportBIN(binFile);
    binFile = main.outputDir + "/optimized/velocity_initial_fluid.bin";
    X0.exportBIN(binFile);
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    std::string binFile = main.outputDir + "/optimized/velocity_initial.bin";
    X0.exportBIN(binFile);
  } else {
    throw std::runtime_error("Undefined gridType");
  }
}
