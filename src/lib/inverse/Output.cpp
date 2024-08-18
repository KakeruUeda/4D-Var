/**
 * @file Output.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "InverseProblem.h"

/***************************************
 * @brief Output velocity and pressure.
 */
void InverseProblem::outputFowardSolutions(const int loop)
{
  if(mpi.myId != 0)
    return;

  if(main.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.updateSolutionsVTI(t);
      main.outputSolutionsVTI("main", t, loop);
    }
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.outputSolutionsVTU("main", t, loop);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }
}

/************************
 * @brief Output w, q, l.
 */
void InverseProblem::outputAdjointSolutions(const int loop)
{
  if(mpi.myId != 0) {
    return;
  }

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

/*******************************************
 * @brief Output control variables X and X0.
 */
void InverseProblem::outputControlVariables(const int loop)
{
  if(mpi.myId != 0) {
    return;
  }

  std::string vtuFile;
  for(int t = 0; t < main.timeMax; t++) {
    vtuFile = main.outputDir + "/other/X_" + to_string(loop) + "_" + to_string(t) + ".vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "X", main.grid.node, main.grid.cell, XArr, t);
  }
  vtuFile = main.outputDir + "/other/X0_" + to_string(loop) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "X", main.grid.node, main.grid.cell, X0Arr);
}

/******************************************
 * @brief Update control variables for VTI.
 */
void InverseProblem::updateControlVariablesVTI()
{
  for(int in = 0; in < main.grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      X0vtiArr(main.grid.vecFluidUniqueNodes[in], d) = X0Arr(in, d);
    }
  }
}

/***********************************************
 * @brief Output velocity data (vCFD, vMRI, ve).
 */
void InverseProblem::outputVelocityData(const int loop)
{
  if(mpi.myId > 0)
    return;

  std::string vtiFile;
  for(int t = 0; t < main.snap.nSnapShot; t++) {
    vtiFile = main.outputDir + "/data/data_" + to_string(loop) + "_" + to_string(t) + ".vti";
    data.exportVTI(vtiFile, t);
  }
}

/*************************************************
 * @brief Output time interpolated feedback force.
 */
void InverseProblem::outputFeedbackForce(const int loop)
{
  if(mpi.myId > 0)
    return;

  std::string vtuFile;
  for(int t = 0; t < main.timeMax; t++) {
    vtuFile = main.outputDir + "/other/feedbackForce" + to_string(loop) + "_" + to_string(t) + ".vtu";
    EXPORT::exportVectorPointDataVTU(vtuFile, "feedbackForce", main.grid.node, main.grid.cell, adjoint.feedbackForceT,
                                     t);
  }
}

/****************************************
 * @brief Output velocity to binary file.
 */
void InverseProblem::outputVelocityBIN(const int loop)
{
  if(mpi.myId > 0)
    return;

  std::string binFile;
  for(int t = 0; t < main.timeMax; t++) {
    main.updateSolutionsVTI(t);
    binFile = main.outputDir + "/bin/velocity_" + to_string(loop) + "_" + to_string(t) + ".bin";
    main.vvti.exportBIN(binFile, t);
  }
  updateControlVariablesVTI();
  binFile = main.outputDir + "/bin/velocity_initial_" + to_string(loop) + ".bin";
  X0vtiArr.exportBIN(binFile);
}

/************************************
 * @brief Output optimized variables.
 */
void InverseProblem::outputOptimizedVariables()
{
  if(mpi.myId > 0)
    return;

  // Optimized velocity field
  if(main.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.updateSolutionsVTI(t);
      main.outputSolutionsVTI("optimized", t);
    }
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.outputSolutionsVTU("optimized", t);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }

  // To binary file
  if(main.grid.gridType == GridType::STRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      main.updateSolutionsVTI(t);
      std::string binFile = main.outputDir + "/optimized/velocity_" + to_string(t) + ".bin";
      main.vvti.exportBIN(binFile, t);
    }
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    for(int t = 0; t < main.timeMax; t++) {
      std::string binFile = main.outputDir + "/optimized/velocity_" + to_string(t) + ".bin";
      main.vt.exportBIN(binFile, t);
    }
  } else {
    throw std::runtime_error("Undefined gridType");
  }

  // Optimized initial velocity field
  if(main.grid.gridType == GridType::STRUCTURED) {
    std::string vtiFile;
    updateControlVariablesVTI();
    vtiFile = main.outputDir + "/optimized/velocity_initial.vti";
    EXPORT::exportVectorPointDataVTI(vtiFile, "velocity_initial", X0vtiArr, main.grid.nx, main.grid.ny, main.grid.nz,
                                     main.grid.dx, main.grid.dy, main.grid.dz);
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    std::string vtiFile;
    vtiFile = main.outputDir + "/optimized/velocity_initial.vtu";
    EXPORT::exportVectorPointDataVTU(vtiFile, "velocity_initial", main.grid.node, main.grid.cell, X0Arr);
  } else {
    throw std::runtime_error("Undefined gridType");
  }

  // TO binary file
  if(main.grid.gridType == GridType::STRUCTURED) {
    updateControlVariablesVTI();
    std::string binFile = main.outputDir + "/optimized/velocity_initial.bin";
    X0vtiArr.exportBIN(binFile);
  } else if(main.grid.gridType == GridType::UNSTRUCTURED) {
    std::string binFile = main.outputDir + "/optimized/velocity_initial.bin";
    X0Arr.exportBIN(binFile);
  } else {
    throw std::runtime_error("Undefined gridType");
  }
}
