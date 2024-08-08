/**
 * @file InverseProblem.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/***************************************************
 * @brief construct inverse object from config param
 */
InverseProblem::InverseProblem(Config &conf):
main(conf), adjoint(conf), data(conf), app(conf.app), 
vvox(conf.vvox), dim(conf.dim), outputDir(conf.outputDir), 
outputItr(conf.outputItr), nOMP(conf.nOMP), aCF(conf.aCF),
bCF(conf.bCF), gCF(conf.gCF), alphaX0(conf.alphaX0), alphaX(conf.alphaX),
loopMax(conf.loopMax), planeDir(conf.planeDir)
{
    std::string dir;
    std::string output = "output";

    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/domain";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/data";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/main";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/adjoint";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/other";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/optimized";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    main.outputDir = outputDir;
    adjoint.outputDir = outputDir;
}

/*****************************
 * @brief Run inverse routine
 */
void InverseProblem::runSimulation()
{ 
    main.outputDomain();
    guessInitialCondition();
    std::ofstream cf(outputDir + "/dat/costFunction.dat");

    for(int loop=0; loop<loopMax; loop++){
        PetscPrintf(MPI_COMM_WORLD, "\nOpt itr. : %d\n", loop);

        compCostFunction();
        PetscPrintf(MPI_COMM_WORLD, "\ncostFunction = %e\n", costFunction.total);

        costFunction.history.push_back(costFunction.total);
        if(mpi.myId == 0){
            cf << costFunction.term1 << " " << costFunction.term2 << " "
               << costFunction.term3 << " " << costFunction.term4 << " "
               << costFunction.term5 << " " << costFunction.term6 << " "
               << costFunction.term7 << " " << costFunction.total;
            cf << std::endl;
        }
        bool isConverged = checkConvergence(cf, loop);
        if(isConverged) break;

        compFeedbackForce();
        compTimeInterpolatedFeedbackForce();
        adjoint.solveAdjoint(main, outputDir, feedbackForceT);

        if(loop % outputItr == 0){
            outputFowardSolutions(loop);
            outputAdjointSolutions(loop);
            outputVelocityBIN(loop);
            outputControlVariables(loop);
            outputVelocityData(loop);
            outputFeedbackForce(loop); 
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        compOptimalCondition();

        armijoCriteriaX0(costFunction.total);
        armijoCriteriaX(costFunction.total);

        PetscPrintf(MPI_COMM_WORLD, "\n(alphaX0, alphaX) = (%f, %f)\n", alphaX0, alphaX);
        updataControlVariables(main);
    }
    cf.close();

    outputOptimizedVariables();
}

/******************************************
 * @brief Check if cost function converged.
 */
bool InverseProblem::checkConvergence(std::ofstream &cf, const int loop)
{    
    if(loop >= 10){
        double diff = costFunction.history[loop] - costFunction.history[loop-10];
        diff = fabs(diff);
        
        if(diff < 1e-4){
            PetscPrintf(MPI_COMM_WORLD, "Converged. OptItr = %d", loop);
            cf.close();
            return true;
        }
    }
    return false;
}

/************************************************************
 * @brief Guess initial condition. 
 *        Getting initial velocity field for inverse problem 
 *        using poiseuille inlet condition.
 */
void InverseProblem::guessInitialCondition()
{
    main.compInitialCondition(main.grid.dirichlet.vDirichletNew, 
                              main.grid.dirichlet.pDirichletNew);

    int snapCount = 0;
    for(int t=0; t<main.timeMax; t++){   
        for(int in=0; in<main.grid.node.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                main.grid.node.vt[t][in][d] = main.grid.node.v0[in][d];
            }
        }
        if((t - main.snap.snapTimeBeginItr) % main.snap.snapInterval == 0){;
            main.snap.takeSnapShot(main.grid.node.vt[t], snapCount, main.grid.node.nNodesGlobal, dim);
            snapCount++;
        }
    }

    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    for(int t=0; t<main.timeMax; t++){
        for(int ib=0; ib<n; ib++){
            int key = adjoint.grid.dirichlet.controlBoundaryMap[ib];
            auto it = main.grid.dirichlet.vDirichlet[t].find(key);
            if(it != main.grid.dirichlet.vDirichlet[t].end()){
                for(int d=0; d<dim; d++){
                    X[t][ib][d] = main.grid.dirichlet.vDirichlet[t][key][d];
                }
            }
        }
    } 
    for(int in=0; in<main.grid.node.nNodesGlobal; in++){
        for(int d=0; d<main.dim; d++){
            X0[in][d] = main.grid.node.v0[in][d];
        }
    }
}

/*****************************************************
 * @brief Update control variables for next iteration.
 */
void InverseProblem::updataControlVariables(DirectProblem &main)
{   
    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    for(int t=0; t<main.timeMax; t++){
        for(int ib=0; ib<n; ib++){
            int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
            std::vector<double> vecTmp(dim, 0e0);
            for(int d=0; d<dim; d++){
                X[t][ib][d] += alphaX * (-grad[t][ib][d]);
            }
        }
    }
    for(int in=0; in<main.grid.node.nNodesGlobal; in++){
        for(int d=0; d<main.dim; d++){
            X0[in][d] += alphaX0 * (-gradInitVel[in][d]);
        }
    }
}
/*****************************************************
 * @brief Update control variables for next iteration.
 */
void InverseProblem::updataControlVariables(DirectProblem &main, const double alphaX, const double alphaX0)
{   
    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    for(int t=0; t<main.timeMax; t++){
        for(int ib=0; ib<n; ib++){
            int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
            std::vector<double> vecTmp(dim, 0e0);
            for(int d=0; d<dim; d++){
                X[t][ib][d] += alphaX * (-grad[t][ib][d]);
            }
        }
    }
    for(int in=0; in<main.grid.node.nNodesGlobal; in++){
        for(int d=0; d<main.dim; d++){
            X0[in][d] += alphaX0 * (-gradInitVel[in][d]);
        }
    }
}

/***************************************
 * @brief Output velocity and pressure.
 */
void InverseProblem::outputFowardSolutions(const int loop)
{
    if(mpi.myId != 0) return;

    for(int t=0; t<main.timeMax; t++){
        switch(main.grid.gridType){
            case GridType::STRUCTURED:
                main.updateSolutionsVTI(t);
                main.outputSolutionsVTI("main", t, loop);
                break;
            case GridType::UNSTRUCTURED:
                main.outputSolutionsVTU("main", t, loop);
            break;
            default:
                PetscPrintf(MPI_COMM_WORLD, "\nUndifined gridType\n"); 
                exit(1);
            break;
        }
    }
    
}

/************************
 * @brief Output w, q, l.
 */
void InverseProblem::outputAdjointSolutions(const int loop)
{
    if(mpi.myId != 0) return;

    for(int t=0; t<adjoint.timeMax; t++){
        switch(adjoint.grid.gridType){
            case GridType::STRUCTURED:
                adjoint.updateSolutionsVTI(t);
                adjoint.outputSolutionsVTI("adjoint", t, loop);
                break;
            case GridType::UNSTRUCTURED:
                adjoint.outputSolutionsVTU("adjoint", t, loop);
            break;
            default:
                PetscPrintf(MPI_COMM_WORLD, "\nUndifined gridType\n"); 
                exit(1);
            break;
        }
    }
}

/*******************************************
 * @brief Output control variables X and X0.
 */
void InverseProblem::outputControlVariables(const int loop)
{
    if(mpi.myId != 0) return;

    std::string vtiFile;
    updateControlVariablesVTI();
    vtiFile = main.outputDir + "/main/X0_" + to_string(loop) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "X0", X0vti, main.grid.nx, main.grid.ny, 
                                  main.grid.nz, main.grid.dx, main.grid.dy, main.grid.dz);
}

/******************************************
 * @brief Update control variables for VTI.
 */
void InverseProblem::updateControlVariablesVTI()
{    
    for(int in=0; in<main.grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            X0vti[main.grid.node.sortNode[in]][d] = X0[in][d];
        }
    }   
}

/***********************************************
 * @brief Output velocity data (vCFD, vMRI, ve).
 */
void InverseProblem::outputVelocityData(const int loop)
{
    if(mpi.myId > 0) return;
    
    std::string vtiFile;
    for(int t=0; t<main.snap.nSnapShot; t++){
       vtiFile = main.outputDir + "/data/data_" + to_string(loop) + "_" + to_string(t) + ".vti";
       VTK::exportVelocityDataVTI(vtiFile, data, t);
    }
}

/*************************************************
 * @brief Output time interpolated feedback force.
 */
void InverseProblem::outputFeedbackForce(const int loop)
{
    if(mpi.myId > 0) return;
    
    std::string vtuFile;
    for(int t=0; t<main.timeMax; t++){
       vtuFile = main.outputDir + "/other/feedbackForce" + to_string(loop) + "_" + to_string(t) + ".vtu";
       VTK::exportVectorPointDataVTU(vtuFile, "feedbackForce", main.grid.node, main.grid.cell, feedbackForceT[t]);
    }
}

/****************************************
 * @brief Output velocity to binary file.
 */
void InverseProblem::outputVelocityBIN(const int loop)
{
    if(mpi.myId > 0) return;

    std::string binFile;
    for(int t=0; t<main.timeMax; t++){
        main.updateSolutionsVTI(t);
        binFile = main.outputDir + "/main/velocity_" + to_string(loop) + "_" + to_string(t) + ".bin";
        BIN::exportVectorDataBIN(binFile, main.grid.node.vvti);
    }
}

/************************************
 * @brief Output optimized variables.
 */
void InverseProblem::outputOptimizedVariables()
{
    if(mpi.myId > 0) return;

    for(int t=0; t<main.timeMax; t++){
        switch(main.grid.gridType){
            case GridType::STRUCTURED:
                main.updateSolutionsVTI(t);
                main.outputSolutionsVTI("optimized", t);
                break;
            case GridType::UNSTRUCTURED:
                main.outputSolutionsVTU("optimized", t);
            break;
            default:
                PetscPrintf(MPI_COMM_WORLD, "\nUndifined gridType\n"); 
                exit(1);
            break;
        }
    }

    std::string vtiFile;
    updateControlVariablesVTI();
    vtiFile = main.outputDir + "/optimized/X0.vti";
    VTK::exportVectorPointDataVTI(vtiFile, "X0", X0vti, main.grid.nx, main.grid.ny, 
                                  main.grid.nz, main.grid.dx, main.grid.dy, main.grid.dz);

    // output bin
    std::string binFile;
    for(int t=0; t<main.timeMax; t++){
        main.updateSolutionsVTI(t);
        binFile = main.outputDir + "/optimized/velocity_" + to_string(t) + ".bin";
        BIN::exportVectorDataBIN(binFile, main.grid.node.vvti);
    }

    binFile = main.outputDir + "/optimized/X0.bin";
    BIN::exportVectorDataBIN(binFile, X0vti);
}

/****************************************************
 * @brief Compute discrepancies over domain and time
 */
void InverseProblem::compCostFunction() 
{  
    costFunction.term1 = 0e0;
    costFunction.term2 = 0e0;
    costFunction.term3 = 0e0;
    costFunction.term4 = 0e0;
    costFunction.term5 = 0e0;
    costFunction.term6 = 0e0;
    costFunction.term7 = 0e0;
    
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<data.nCellsGlobal; ic++){
            for(int d=0; d<dim; d++){           
                data(ic).vCFD[t][d] = 0e0;
                data(ic).ve[t][d] = 0e0;
            }
        }
    }

    switch(vvox){
        case VoxelVelocity::AVERAGE:
            for(int t=0; t<main.snap.nSnapShot; t++){
                for(int ic=0; ic<data.nCellsGlobal; ic++){
                    data(ic).average(main.grid.cell, main.snap.v[t], t, main.dim);
                }
            }
        break;

        case VoxelVelocity::INTERPOLATION:
            for(int t=0; t<main.snap.nSnapShot; t++){
                for(int ic=0; ic<data.nCellsGlobal; ic++){
                    data(ic).interpolate(main.grid.node, main.grid.cell, main.snap.v[t], t, main.dim);
                }
            }
        break;
    
        default:
            PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
            exit(1);
    }

    // term 1
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<data.nCellsGlobal; ic++){
            double dev = 0e0;
            double volume = data.dx * data.dy * data.dz;
            double deltaT = main.dt * main.snap.snapInterval;
            for(int d=0; d<dim; d++){
                data(ic).ve[t][d] = data(ic).vCFD[t][d] - data(ic).vMRI[t][d];
                dev += data(ic).ve[t][d] * data(ic).ve[t][d];
            }
            costFunction.term1 += 5e-1 * aCF * dev * volume * deltaT;
        }
    }

    int nc = adjoint.grid.dirichlet.nControlNodesInCell;

    Gauss gauss(2);
    Function func2d(nc, dim-1);

    // term 2, 3, 4, 5
    for(int t=0; t<main.timeMax; t++){
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            for(int p=0; p<nc; p++){
                for(int d=0; d<dim-1; d++){
                    int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    func2d.xCurrent[p][d] = main.grid.node.x[n][planeDir[d]];
                }
            }

            double value2 = 0e0;
            double value3 = 0e0;
            double value4 = 0e0;
            double value5 = 0e0;

            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    func2d.weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(func2d.N, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(func2d.dNdr, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralRegTerm2(func2d, value2, ic, t);
                    GaussIntegralRegTerm3(func2d, value3, ic, t);
                    GaussIntegralRegTerm4(func2d, value4, ic, t);
                    GaussIntegralRegTerm5(func2d, value5, ic, t);
                }
            }
            costFunction.term2 += 5e-1 * bCF * value2;
            costFunction.term3 += 5e-1 * bCF * value3;
            costFunction.term4 += 5e-1 * bCF * value4;
            costFunction.term5 += 5e-1 * bCF * value5;
        }
    }

    Function func3d(main.grid.cell.nNodesInCell, dim);

    // term 6, 7
    for(int ic=0; ic<main.grid.cell.nCellsGlobal; ic++){
        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            for(int d=0; d<main.dim; d++){
                func3d.xCurrent[p][d] = main.grid.node.x[main.grid.cell(ic).node[p]][d];
            }
        }

        double value6 = 0e0;
        double value7 = 0e0;

        for(int i1=0; i1<2; i1++){
            for(int i2=0; i2<2; i2++){
                for(int i3=0; i3<2; i3++){
                    func3d.weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                    ShapeFunction3D::C3D8_N(func3d.N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    ShapeFunction3D::C3D8_dNdr(func3d.dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    GaussIntegralRegTerm6(func3d, value6, ic);
                    GaussIntegralRegTerm7(func3d, value7, ic);
                }
            }
        }
        costFunction.term6 += 5e-1 * gCF * value6;
        costFunction.term7 += 5e-1 * gCF * value7;
    }

    costFunction.sum();
}

/************************************************
 * @brief Compute value for regularization term2
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralRegTerm2(Function &func, double &value, const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 

    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<double> u(dim, 0e0);

    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<dim; d++){
            u[d] += func.N[p] * main.grid.node.vt[t][n][d];
        }
    }

    for(int d=0; d<dim; d++){
        value += u[d] * u[d] * func.vol;
    }
}

/************************************************
 * @brief Compute value for regularization term3 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralRegTerm3(Function &func, double &value, const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 
    
    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    MathFEM::comp_dNdx2D(func.dNdx, func.dNdr, dxdr, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<std::vector<double>> dudx;
    VecTool::resize(dudx, dim, dim-1);

    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                dudx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t][n][d1];
            }
        }
    }

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            value += dudx[d1][d2] * dudx[d1][d2] * func.vol;
        }
    }
}

/************************************************
 * @brief Compute value for regularization term4
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralRegTerm4(Function &func, double &value, const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 

    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<double> u(dim, 0e0);
    std::vector<double> ub(dim, 0e0);
    std::vector<double> dudt(dim, 0e0);

    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<dim; d++){
            if(t == 0){
                u[d] += func.N[p] * main.grid.node.vt[t][n][d];
                ub[d] += func.N[p] * main.grid.node.v0[n][d];
            }else{
                u[d] += func.N[p] * main.grid.node.vt[t][n][d];
                ub[d] += func.N[p] * main.grid.node.vt[t-1][n][d];
            }
        }
    }

    for(int d=0; d<dim; d++){
        dudt[d] = (u[d] - ub[d]) / main.dt;
    }
    for(int d=0; d<dim; d++){
        value += dudt[d] * dudt[d] * func.vol;
    }
}

/************************************************
 * @brief Compute value for regularization term3 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralRegTerm5(Function &func, double &value, const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 
    
    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    MathFEM::comp_dNdx2D(func.dNdx, func.dNdr, dxdr, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<std::vector<double>> dudx;
    std::vector<std::vector<double>> dubdx;
    std::vector<std::vector<double>> dudxdt;

    VecTool::resize(dudx, dim, dim-1);
    VecTool::resize(dubdx, dim, dim-1);
    VecTool::resize(dudxdt, dim, dim-1);
    
    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                if(t == 0){
                    dudx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t][n][d1];
                    dubdx[d1][d2] += func.dNdx[p][d2] * main.grid.node.v0[n][d1];
                }else{
                    dudx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t][n][d1];
                    dubdx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t-1][n][d1];
                }
            }
        }
    }

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            dudxdt[d1][d2] = (dudx[d1][d2] - dubdx[d1][d2]) / main.dt;
        }
    }
    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            value += dudxdt[d1][d2] * dudxdt[d1][d2] * func.vol;
        }
    }
}

/************************************************
 * @brief Compute value for regularization term6
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralRegTerm6(Function &func, double &value, const int ic)
{
    double dxdr[3][3];

    MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, main.grid.cell.nNodesInCell);
    func.detJ = MathCommon::compDeterminant_3x3(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<double> u0(dim, 0e0);

    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        int n = main.grid.cell(ic).node[p];
        for(int d=0; d<dim; d++){
            u0[d] += func.N[p] * main.grid.node.v0[n][d];
        }
    }

    for(int d=0; d<dim; d++){
        value += u0[d] * u0[d] * func.vol;
    }
}

/************************************************
 * @brief Compute value for regularization term7
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralRegTerm7(Function &func, double &value, const int ic)
{
    double dxdr[3][3];

    MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, main.grid.cell.nNodesInCell);
    MathFEM::comp_dNdx(func.dNdx, func.dNdr, dxdr, main.grid.cell.nNodesInCell);
    func.detJ = MathCommon::compDeterminant_3x3(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<std::vector<double>> du0dx;
    VecTool::resize(du0dx, dim, dim);
    
    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        int n = main.grid.cell(ic).node[p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<3; d2++){
                du0dx[d1][d2] += func.dNdx[p][d2] * main.grid.node.v0[n][d1];
            }
        }
    }

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<3; d2++){
            value += du0dx[d1][d2] * du0dx[d1][d2] * func.vol;
        }
    }
}

/************************************************
 * @brief Compute value for regularization term3 
 *        on gauss integral point.
 */
void InverseProblem::compFeedbackForce()
{   
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int in=0; in<main.grid.node.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                feedbackForce[t][in][d] = 0e0;
            }
        }
        data.compEdgeValue(t);
        for(int ic=0; ic<main.grid.cell.nCellsGlobal; ic++){
            Function func3d(adjoint.grid.cell.nNodesInCell, dim);
            assembleFeedbackForce(func3d, ic, t);
        }
    }
}

/********************************************************
 * @brief Assemble RHS feedback force for adjoint system.
 */
void InverseProblem::assembleFeedbackForce(Function &func, const int ic, const int t)
{
    double dxdr[3][3];
    Gauss g2(2); 

    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        for(int d=0; d<dim; d++){
            func.xCurrent[p][d] = main.grid.cell(ic).x[p][d];
        }
    }
    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            for(int i3=0; i3<2; i3++){
                ShapeFunction3D::C3D8_N(func.N, g2.point[i1], g2.point[i2], g2.point[i3]);
                ShapeFunction3D::C3D8_dNdr(func.dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);
                MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, main.grid.cell.nNodesInCell);

                func.detJ = MathCommon::compDeterminant_3x3(dxdr);
                func.weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];
                double feedback[3] = {0e0, 0e0, 0e0};
                double point[3] = {0e0, 0e0, 0e0};

                for(int d=0; d<dim; d++){
                    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
                        point[d] += func.N[p] * func.xCurrent[p][d];
                    }
                }
                compInterpolatedFeeback(feedback, point);
                feedbackGaussIntegral(func, feedback, ic, t);
            }
        }
    }
}

/*****************************************************************
 * @brief Interpolate space-discrete feedback force onto CFD node.
 */
void InverseProblem::compInterpolatedFeeback(double (&feedback)[3], double (&point)[3])
{
    double px = point[0] + (data.dx / 2e0);
    double py = point[1] + (data.dy / 2e0);
    double pz = point[2] + (data.dz / 2e0);

    int ix = (px / data.dx) + EPS;
    int iy = (py / data.dy) + EPS;
    int iz = (pz / data.dz) + EPS;

    double s = (px - (ix * data.dx + (data.dx / 2e0)));
    double t = (py - (iy * data.dy + (data.dy / 2e0)));
    double u = (pz - (iz * data.dz + (data.dz / 2e0)));

    s = s / (data.dx / 2e0);
    t = t / (data.dy / 2e0);
    u = u / (data.dz / 2e0);

    if(s<-1-EPS || s>1+EPS){
        PetscPrintf(MPI_COMM_WORLD, "\ns interpolation error found.\n");
    }else if(t<-1-EPS || t>1+EPS){
        PetscPrintf(MPI_COMM_WORLD, "\nt interpolation error found.\n");
    }else if(u<-1-EPS || u>1+EPS){
        PetscPrintf(MPI_COMM_WORLD, "\nu interpolation error found.\n");
    }

    std::vector<double> N;
    VecTool::resize(N, main.grid.cell.nNodesInCell);
    ShapeFunction3D::C3D8_N(N, s, t, u);

    for(int d=0; d<dim; d++){
        feedback[d] = N[0]*data.vEX[iz][iy][ix][d]       +  N[1]*data.vEX[iz][iy][ix+1][d]
                    + N[2]*data.vEX[iz][iy+1][ix+1][d]   +  N[3]*data.vEX[iz][iy+1][ix][d]
                    + N[4]*data.vEX[iz+1][iy][ix][d]     +  N[5]*data.vEX[iz+1][iy][ix+1][d]
                    + N[6]*data.vEX[iz+1][iy+1][ix+1][d] +  N[7]*data.vEX[iz+1][iy+1][ix][d];
    } 
}

/********************************************************
 * @brief Compute feedback force on gauss integral point.
 */
void InverseProblem::feedbackGaussIntegral(Function &func, double (&feedback)[3], const int ic, const int t)
{
    func.vol = func.detJ * func.weight;
    for(int d=0; d<dim; d++){
        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            int in = main.grid.cell(ic).node[p];
            feedbackForce[t][in][d] += aCF * func.N[p] * feedback[d] * func.vol;
        }
    }
}

/**************************************************************
 * @brief Interpolate the feedback force at each CFD time step.
 */
void InverseProblem::compTimeInterpolatedFeedbackForce()
{
    for(int in=0; in<adjoint.grid.node.nNodesGlobal; in++){
        std::vector<double> x, y1, y2, y3;

        for(int t=0; t<main.snap.nSnapShot; t++){
            double dp = t * main.dt * main.snap.snapInterval;
            x.push_back(dp);
            y1.push_back(feedbackForce[t][in][0]);
            y2.push_back(feedbackForce[t][in][1]);  
            y3.push_back(feedbackForce[t][in][2]); 
        }

        vector<Spline::Coefficients> cf_x = Spline::compCoefficients(x, y1);
        vector<Spline::Coefficients> cf_y = Spline::compCoefficients(x, y2);
        vector<Spline::Coefficients> cf_z = Spline::compCoefficients(x, y3);
        
        for(int t=0; t<adjoint.timeMax; t++){
            double p = t * main.dt;
            feedbackForceT[t][in][0] = Spline::evaluate(cf_x, p);
            feedbackForceT[t][in][1] = Spline::evaluate(cf_y, p);
            feedbackForceT[t][in][2] = Spline::evaluate(cf_z, p);
        }
    }

}

/**********************************************
 * @brief Compute gradient of control variables.
 */
void InverseProblem::compOptimalCondition()
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    Function func2d(nc, dim-1);
    Gauss gauss(2);

    // Optimal condition for inlet boundary condition
    for(int t=0; t<main.timeMax; t++){
        for(int in=0; in<adjoint.grid.node.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                gradWholeNode[t][in][d] = 0e0;
            }
        }        
        for(int ib=0; ib<adjoint.grid.dirichlet.controlBoundaryMap.size(); ib++){
            for(int d=0; d<dim; d++){
                grad[t][ib][d] = 0e0;
            }
        }

        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            std::vector<std::vector<double>> value1, value2, value3, value4, value5;
        
            VecTool::resize(value1, main.grid.cell.nNodesInCell, dim);
            VecTool::resize(value2, main.grid.cell.nNodesInCell, dim);
            VecTool::resize(value3, main.grid.cell.nNodesInCell, dim);
            VecTool::resize(value4, main.grid.cell.nNodesInCell, dim);
            VecTool::resize(value5, main.grid.cell.nNodesInCell, dim);

            for(int p=0; p<nc; p++){
                for(int d=0; d<dim-1; d++){
                    int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    func2d.xCurrent[p][d] = main.grid.node.x[in][planeDir[d]];
                }
            }
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    func2d.weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(func2d.N, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(func2d.dNdr, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralOptimalConditionXTerm1(func2d, value1, ic, t);
                    GaussIntegralOptimalConditionXTerm2(func2d, value2, ic, t);
                    GaussIntegralOptimalConditionXTerm3(func2d, value3, ic, t);
                    GaussIntegralOptimalConditionXTerm4(func2d, value4, ic, t);
                    GaussIntegralOptimalConditionXTerm5(func2d, value5, ic, t);
                }
            }
            for(int p=0; p<nc; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += bCF * value1[p][d];
                    gradWholeNode[t][in][d] += bCF * value2[p][d];
                    gradWholeNode[t][in][d] += bCF * value3[p][d];
                    gradWholeNode[t][in][d] += bCF * value4[p][d];
                    gradWholeNode[t][in][d] +=       value5[p][d];
                }
            }
        }
        
        for(int ib=0; ib<adjoint.grid.dirichlet.controlBoundaryMap.size(); ib++){
            int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
            for(int d=0; d<dim; d++){
                grad[t][ib][d] = gradWholeNode[t][in][d];
            }
        }
    }

    // Optimal condition for initial velocity field
    Function func3d(main.grid.cell.nNodesInCell, dim);

    for(int in=0; in<main.grid.node.nNodesGlobal; in++){
        for(int d=0; d<main.dim; d++){
            gradInitVel[in][d] = 0e0;
        }
    }

    for(int ic=0; ic<main.grid.cell.nCellsGlobal; ic++){
        std::vector<std::vector<double>> value1, value2, value3;

        VecTool::resize(value1, main.grid.cell.nNodesInCell, dim);
        VecTool::resize(value2, main.grid.cell.nNodesInCell, dim);
        VecTool::resize(value3, main.grid.cell.nNodesInCell, dim);

        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            for(int d=0; d<main.dim; d++){
                func3d.xCurrent[p][d] = main.grid.node.x[main.grid.cell(ic).node[p]][d];
            }
        }
        for(int i1=0; i1<2; i1++){
            for(int i2=0; i2<2; i2++){
                for(int i3=0; i3<2; i3++){
                    func3d.weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                    ShapeFunction3D::C3D8_N(func3d.N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    ShapeFunction3D::C3D8_dNdr(func3d.dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    GaussIntegralOptimalConditionX0Term1(func3d, value1, ic);
                    GaussIntegralOptimalConditionX0Term2(func3d, value2, ic);
                    GaussIntegralOptimalConditionX0Term3(func3d, value3, ic);
                }
            }
        }

        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            int in = main.grid.cell(ic).node[p];
            for(int d=0; d<main.dim; d++){
                gradInitVel[in][d] += gCF * value1[p][d];
                gradInitVel[in][d] += gCF * value2[p][d];
                gradInitVel[in][d] +=       value3[p][d];
            }
        }
    }
}

/***********************************************
 * @brief Set values needed for matrix assembly 
 *        on gauss integral points.
 */
void InverseProblem::setValue(Function &func, const int ic)
{
    // main var - v
    for(int d=0; d<main.dim; d++){
        adjoint.vk[d] = 0e0;
        adjoint.vk1[d] = 0e0;
        adjoint.vk2[d] = 0e0;
        for(int p=0; p<adjoint.grid.cell.nNodesInCell; p++){
            int n = adjoint.grid.cell(ic).node[p];
            adjoint.vk[d]  += func.N[p] * main.grid.node.vt[0][n][d];
            adjoint.vk1[d] += func.N[p] * main.grid.node.vt[1][n][d];
            adjoint.vk2[d] += func.N[p] * main.grid.node.vt[2][n][d];
        }
    }

    // main var - dvdx
    for(int d=0; d<main.dim; d++){
        for(int e=0; e<main.dim; e++){
            adjoint.dvkdx[d][e] = 0e0;
            adjoint.dvk1dx[d][e] = 0e0;
            adjoint.dvk2dx[d][e] = 0e0;
            for(int p=0; p<adjoint.grid.cell.nNodesInCell; p++){
                int n = adjoint.grid.cell(ic).node[p];
                adjoint.dvkdx[d][e]  += func.dNdx[p][e] * main.grid.node.vt[0][n][d];
                adjoint.dvk1dx[d][e] += func.dNdx[p][e] * main.grid.node.vt[1][n][d];
                adjoint.dvk2dx[d][e] += func.dNdx[p][e] * main.grid.node.vt[2][n][d];
            }
        }
    }

    // main var - dpdx
    for(int d=0; d<main.dim; d++){
        adjoint.dpkdx[d] = 0e0;
        adjoint.dpk1dx[d] = 0e0;
        adjoint.dpk2dx[d] = 0e0;
        for(int p=0; p<adjoint.grid.cell.nNodesInCell; p++){
            int n = adjoint.grid.cell(ic).node[p];
            adjoint.dpkdx[d]  += func.dNdx[p][d] * main.grid.node.pt[0][n];
            adjoint.dpk1dx[d] += func.dNdx[p][d] * main.grid.node.pt[1][n];
            adjoint.dpk2dx[d] += func.dNdx[p][d] * main.grid.node.pt[2][n];
        }
    }

    // main var - adv
    for(int d=0; d<main.dim; d++){
        adjoint.advk1[d] = 0e0; 
        adjoint.advk2[d] = 0e0;
        adjoint.advk3[d] = 0e0;
        for(int p=0; p<adjoint.grid.cell.nNodesInCell; p++){
            int n = adjoint.grid.cell(ic).node[p];
            adjoint.advk1[d] += func.N[p] * main.grid.node.v0[n][d];
            adjoint.advk2[d] += func. N[p] * (1.5 * main.grid.node.vt[0][n][d]
                              - 0.5 * main.grid.node.v0[n][d]);
            adjoint.advk3[d] += func.N[p] * (1.5 * main.grid.node.vt[1][n][d] 
                              - 0.5 * main.grid.node.vt[0][n][d]);
        }
    }

    // lagrange multiplier - w
    for(int d=0; d<main.dim; d++){
        adjoint.wk1[d] = 0e0;
        adjoint.wk2[d] = 0e0;
        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            int n = adjoint.grid.cell(ic).node[p];
            adjoint.wk1[d] += func.N[p] * adjoint.grid.node.wt[0][n][d];
            adjoint.wk2[d] += func.N[p] * adjoint.grid.node.wt[1][n][d];
        }
    }

    // lagrange multiplier - dwdx
    for(int d=0; d<main.dim; d++){
        for(int e=0; e<main.dim; e++){
            adjoint.dwk1dx[d][e] = 0e0;
            adjoint.dwk2dx[d][e] = 0e0;
            for(int p=0; p<main.grid.cell.nNodesInCell; p++){
                int n = adjoint.grid.cell(ic).node[p];
                adjoint.dwk1dx[d][e] += func.dNdx[p][e] * adjoint.grid.node.wt[0][n][d];
                adjoint.dwk2dx[d][e] += func.dNdx[p][e] * adjoint.grid.node.wt[1][n][d];
            }
        }
    }

    // lagrange multiplier - dqdx
    for(int d=0; d<main.dim; d++){
        adjoint.dqk1dx[d] = 0e0;
        adjoint.dqk2dx[d] = 0e0;
        for(int p=0; p<adjoint.grid.cell.nNodesInCell; p++){
            int n = adjoint.grid.cell(ic).node[p];
            adjoint.dqk1dx[d] += func.dNdx[p][d] * adjoint.grid.node.qt[0][n];
            adjoint.dqk2dx[d] += func.dNdx[p][d] * adjoint.grid.node.qt[1][n];
        }
    }

}

/*****************************************************
 * @brief Compute value for term1 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionXTerm1(Function &func, std::vector<std::vector<double>> &value, 
                                                         const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 
    
    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    double vel[3] = {0e0, 0e0, 0e0};

    for(int d=0; d<3; d++){
        vel[d] = 0e0;
    }
    for(int p=0; p<nc; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<3; d++){
            vel[d] += func.N[p] * main.grid.node.vt[t][in][d];
        }
    }
    for(int p=0; p<nc; p++){
        for(int d=0; d<3; d++){
            value[p][d] += vel[d] * func.N[p] * func.vol;
        }
    }
}

/*****************************************************
 * @brief Compute value for term2 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionXTerm2(Function &func, std::vector<std::vector<double>> &value, 
                                                         const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 

    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    MathFEM::comp_dNdx2D(func.dNdx, func.dNdr, dxdr, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    double dudx[3][2];

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            dudx[d1][d2] = 0e0;
        }
    } 
    for(int p=0; p<nc; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                dudx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t][in][d1];
            }
        }
    }
    for(int p=0; p<nc; p++){
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                value[p][d1] += dudx[d1][d2] * func.N[p]* func.vol;
            }
        }
    }
}

/*****************************************************
 * @brief Compute value for term3 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionXTerm3(Function &func, std::vector<std::vector<double>> &value, 
                                                         const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 

    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<double> u(dim, 0e0);
    std::vector<double> ub(dim, 0e0);
    std::vector<double> dudt(dim, 0e0);

    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<dim; d++){
            if(t == 0){
                u[d] += func.N[p] * main.grid.node.vt[t][n][d];
                ub[d] += func.N[p] * main.grid.node.v0[n][d];
            }else{
                u[d] += func.N[p] * main.grid.node.vt[t][n][d];
                ub[d] += func.N[p] * main.grid.node.vt[t-1][n][d];
            }
        }
    }

    for(int d=0; d<dim; d++){
        dudt[d] = (u[d] - ub[d]) / main.dt;
    }

    for(int p=0; p<nc; p++){
        for(int d=0; d<3; d++){
            value[p][d] += dudt[d] * func.N[p] * func.vol;
        }
    }
}

/*****************************************************
 * @brief Compute value for term4 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionXTerm4(Function &func, std::vector<std::vector<double>> &value, 
                                                         const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 
    
    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    MathFEM::comp_dNdx2D(func.dNdx, func.dNdr, dxdr, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<std::vector<double>> dudx;
    std::vector<std::vector<double>> dubdx;
    std::vector<std::vector<double>> dudxdt;

    VecTool::resize(dudx, dim, dim-1);
    VecTool::resize(dubdx, dim, dim-1);
    VecTool::resize(dudxdt, dim, dim-1);
    
    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                if(t == 0){
                    dudx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t][n][d1];
                    dubdx[d1][d2] += func.dNdx[p][d2] * main.grid.node.v0[n][d1];
                }else{
                    dudx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t][n][d1];
                    dubdx[d1][d2] += func.dNdx[p][d2] * main.grid.node.vt[t-1][n][d1];
                }
            }
        }
    }

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            dudxdt[d1][d2] = (dudx[d1][d2] - dubdx[d1][d2]) / main.dt;
        }
    }

    for(int p=0; p<nc; p++){
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                value[p][d1] += dudxdt[d1][d2] * func.N[p] * func.vol;
            }
        }
    }
}

/*****************************************************
 * @brief Compute value for term5 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionXTerm5(Function &func, std::vector<std::vector<double>> &value, 
                                                         const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 

    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    double lgp[3];
    
    for(int d=0; d<3; d++){
        lgp[d] = 0e0;
    }
    for(int p=0; p<nc; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<3; d++){
            lgp[d] -= func.N[p] * adjoint.grid.node.lt[t][n][d];
        }
    }
    for(int p=0; p<nc; p++){
        for(int d=0; d<3; d++){
            value[p][d] += lgp[d] * func.N[p] * func.vol;
        }
    }
}

/*****************************************************
 * @brief Compute value for term1 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionX0Term1(Function &func, std::vector<std::vector<double>> &value, const int ic)
{
    double dxdr[3][3];

    MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, main.grid.cell.nNodesInCell);
    func.detJ = MathCommon::compDeterminant_3x3(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<double> u0(dim, 0e0);

    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        int n = main.grid.cell(ic).node[p];
        for(int d=0; d<dim; d++){
            u0[d] += func.N[p] * main.grid.node.v0[n][d];
        }
    }

    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] += u0[d] * func.N[p] * func.vol;
        }
    }
}

/*****************************************************
 * @brief Compute value for term2 in optimal condition 
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionX0Term2(Function &func, std::vector<std::vector<double>> &value, const int ic)
{
    double dxdr[3][3];

    MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, main.grid.cell.nNodesInCell);
    MathFEM::comp_dNdx(func.dNdx, func.dNdr, dxdr, main.grid.cell.nNodesInCell);
    func.detJ = MathCommon::compDeterminant_3x3(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<std::vector<double>> du0dx;
    VecTool::resize(du0dx, dim, dim);
    
    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        int n = main.grid.cell(ic).node[p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<3; d2++){
                du0dx[d1][d2] += func.dNdx[p][d2] * main.grid.node.v0[n][d1];
            }
        }
    }

    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<3; d2++){
                value[p][d1] += du0dx[d1][d2] * func.N[p]* func.vol;
            }
        }
    }
}

/*****************************************************
 * @brief Compute value for term3 in optimal condition
 *        on gauss integral point.
 */
void InverseProblem::GaussIntegralOptimalConditionX0Term3(Function &func, std::vector<std::vector<double>> &value, const int ic)
{
    setValue(func, ic);
    double he = fabs(func.xCurrent[1][0] - func.xCurrent[0][0]);

    adjoint.tau = MathFEM::comp_tau(adjoint.advk2, he, main.Re, main.dt);
    double dxdr[3][3];

    MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, main.grid.cell.nNodesInCell);
    MathFEM::comp_dNdx(func.dNdx, func.dNdr, dxdr, main.grid.cell.nNodesInCell);
    
    func.detJ = MathCommon::compDeterminant_3x3(dxdr);
    func.vol = func.detJ * func.weight;

    int n1, n2, n3;
    func.vol = func.detJ * func.weight;
    double f = main.resistance * main.alpha * (1e0 - main.grid.cell(ic).phi) 
                               / (main.alpha + main.grid.cell(ic).phi);

    for(int p=0; p<main.grid.cell.nNodesInCell; p++){
        // Mass term
        value[p][0] -= func.N[p] * adjoint.wk1[0] / main.dt * func.vol;
        value[p][1] -= func.N[p] * adjoint.wk1[1] / main.dt * func.vol;
        value[p][2] -= func.N[p] * adjoint.wk1[2] / main.dt * func.vol;

        /*
        // Diffusion term
        for(int d=0; d<3; d++){
            if(d == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
            if(d == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
            if(d == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
            value[p][0] += 5e-1 * n1 * func.dNdx[p][d] * adjoint.dwk1dx[0][d] / main.Re * func.vol;
            value[p][1] += 5e-1 * n2 * func.dNdx[p][d] * adjoint.dwk1dx[1][d] / main.Re * func.vol;
            value[p][2] += 5e-1 * n3 * func.dNdx[p][d] * adjoint.dwk1dx[2][d] / main.Re * func.vol;
        }
        value[p][0] += 5e-1 * func.dNdx[p][1] * adjoint.dwk1dx[1][0] / main.Re * func.vol;
        value[p][0] += 5e-1 * func.dNdx[p][2] * adjoint.dwk1dx[2][0] / main.Re * func.vol;
        value[p][1] += 5e-1 * func.dNdx[p][0] * adjoint.dwk1dx[0][1] / main.Re * func.vol;
        value[p][1] += 5e-1 * func.dNdx[p][2] * adjoint.dwk1dx[2][1] / main.Re * func.vol;
        value[p][2] += 5e-1 * func.dNdx[p][0] * adjoint.dwk1dx[0][2] / main.Re * func.vol;
        value[p][2] += 5e-1 * func.dNdx[p][1] * adjoint.dwk1dx[1][2] / main.Re * func.vol;
        */
 
        // Diffusion term
        for(int d=0; d<3; d++){
            value[p][0] += 5e-1 * func.dNdx[p][d] * adjoint.dwk1dx[0][d] / main.Re * func.vol;
            value[p][1] += 5e-1 * func.dNdx[p][d] * adjoint.dwk1dx[1][d] / main.Re * func.vol;
            value[p][2] += 5e-1 * func.dNdx[p][d] * adjoint.dwk1dx[2][d] / main.Re * func.vol;
        }

        // Advection term
        value[p][0] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][0] * adjoint.wk1[0] * func.vol;
        value[p][0] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][0] * adjoint.wk1[0] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][0] += 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * adjoint.wk1[0] * func.vol;
        }
        value[p][0] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][0] * adjoint.wk1[1] * func.vol;
        value[p][0] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][0] * adjoint.wk1[1] * func.vol;
        value[p][0] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][0] * adjoint.wk1[2] * func.vol;
        value[p][0] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][0] * adjoint.wk1[2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][0] += 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][0] * adjoint.wk2[d] * func.vol;
            value[p][0] += 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][0] * adjoint.wk2[d] * func.vol;
        }

        value[p][1] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][1] * adjoint.wk1[1] * func.vol;
        value[p][1] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][1] * adjoint.wk1[1] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][1] += 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * adjoint.wk1[1] * func.vol;
        }
        value[p][1] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][1] * adjoint.wk1[0] * func.vol;
        value[p][1] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][1] * adjoint.wk1[0] * func.vol;
        value[p][1] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][1] * adjoint.wk1[2] * func.vol;
        value[p][1] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][1] * adjoint.wk1[2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][1] += 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][1] * adjoint.wk2[d] * func.vol;
            value[p][1] += 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][1] * adjoint.wk2[d] * func.vol;
        }

        value[p][2] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][2] * adjoint.wk1[2] * func.vol;
        value[p][2] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][2] * adjoint.wk1[2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][2] += 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * adjoint.wk1[2] * func.vol;
        }
        value[p][2] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][2] * adjoint.wk1[0] * func.vol;
        value[p][2] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][2] * adjoint.wk1[0] * func.vol;
        value[p][2] += 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][2] * adjoint.wk1[1] * func.vol;
        value[p][2] += 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][2] * adjoint.wk1[1] * func.vol;   
        for(int d=0; d<main.dim; d++){
            value[p][2] += 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][2] * adjoint.wk2[d] * func.vol;
            value[p][2] += 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][2] * adjoint.wk2[d] * func.vol;
        }

        // Darcy term
        value[p][0] += 5e-1 * f * func.N[p] * adjoint.wk1[0] * func.vol;
        value[p][1] += 5e-1 * f * func.N[p] * adjoint.wk1[1] * func.vol;
        value[p][2] += 5e-1 * f * func.N[p] * adjoint.wk1[2] * func.vol;

        // SUPG mass term
        for(int d=0; d<main.dim; d++){
            value[p][0] += adjoint.tau * 1.5 * func.N[p] * (adjoint.vk1[d] - adjoint.vk[d]) / main.dt * adjoint.dwk1dx[d][0] * func.vol;
            value[p][0] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dwk2dx[d][0] * (adjoint.vk2[d] - adjoint.vk1[d]) / main.dt * func.vol;
        }
        value[p][0] -= adjoint.tau * func.N[p] * adjoint.advk2[0] * adjoint.dwk1dx[0][0] / main.dt * func.vol;

        for(int d=0; d<main.dim; d++){
            value[p][1] += adjoint.tau * 1.5 * func.N[p] * (adjoint.vk1[d] - adjoint.vk[d]) / main.dt * adjoint.dwk1dx[d][1] * func.vol;
            value[p][1] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dwk2dx[d][1] * (adjoint.vk2[d] - adjoint.vk1[d]) / main.dt * func.vol;
        }
        value[p][1] -= adjoint.tau * func.N[p] * adjoint.advk2[1] * adjoint.dwk1dx[1][1] / main.dt * func.vol;

        for(int d=0; d<main.dim; d++){
            value[p][2] += adjoint.tau * 1.5 * func.N[p] * (adjoint.vk1[d] - adjoint.vk[d]) / main.dt * adjoint.dwk1dx[d][1] * func.vol;
            value[p][2] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dwk2dx[d][2] * (adjoint.vk2[d] - adjoint.vk1[d]) / main.dt * func.vol;
        }
        value[p][2] -= adjoint.tau * func.N[p] * adjoint.advk2[2] * adjoint.dwk1dx[2][2] / main.dt * func.vol;

        // SUPG advection term
        std::vector<double> frontAdv2, frontAdv3;
        VecTool::resize(frontAdv2, main.dim);
        VecTool::resize(frontAdv3, main.dim);

        for(int d1=0; d1<main.dim; d1++){
            for(int d2=0; d2<main.dim; d2++){
                frontAdv2[d1] += adjoint.advk2[d2] * adjoint.dwk1dx[d1][d2];
                frontAdv3[d1] += adjoint.advk3[d2] * adjoint.dwk2dx[d1][d2];
            }
        }

        value[p][0] += adjoint.tau * frontAdv2[0] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][0] * func.vol;
        value[p][0] += adjoint.tau * frontAdv2[0] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][0] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][0] += adjoint.tau * frontAdv2[0] * 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * func.vol;
        }
        value[p][0] += adjoint.tau * frontAdv2[1] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][0] * func.vol;
        value[p][0] += adjoint.tau * frontAdv2[1] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][0] * func.vol;
        value[p][0] += adjoint.tau * frontAdv2[2] *  0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][0] * func.vol;
        value[p][0] += adjoint.tau * frontAdv2[2] *  0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][0] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][0] += adjoint.tau * frontAdv3[d] * 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][0] * func.vol;
            value[p][0] += adjoint.tau * frontAdv3[d] * 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][0] * func.vol;
        }

        value[p][1] += adjoint.tau * frontAdv2[1] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][1] * func.vol;
        value[p][1] += adjoint.tau * frontAdv2[1] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][1] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][1] += adjoint.tau * frontAdv2[1] * 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * func.vol;
        }
        value[p][1] += adjoint.tau * frontAdv2[0] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][1] * func.vol;
        value[p][1] += adjoint.tau * frontAdv2[0] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][1] * func.vol;
        value[p][1] += adjoint.tau * frontAdv2[2] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][1] * func.vol;
        value[p][1] += adjoint.tau * frontAdv2[2] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][1] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][1] += adjoint.tau * frontAdv3[d] * 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][1] * func.vol;
            value[p][1] += adjoint.tau * frontAdv3[d] * 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][1] * func.vol;
        }

        value[p][2] += adjoint.tau * frontAdv2[2] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][2] * func.vol;
        value[p][2] += adjoint.tau * frontAdv2[2] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][2] += adjoint.tau * frontAdv2[2] * 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * func.vol;
        }
        value[p][2] += adjoint.tau * frontAdv2[0] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][2] * func.vol;
        value[p][2] += adjoint.tau * frontAdv2[0] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][2] * func.vol;
        value[p][2] += adjoint.tau * frontAdv2[1] * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][2] * func.vol;
        value[p][2] += adjoint.tau * frontAdv2[1] * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][2] * func.vol;   
        for(int d=0; d<main.dim; d++){
            value[p][2] += adjoint.tau * frontAdv3[d] * 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][2] * func.vol;
            value[p][2] += adjoint.tau * frontAdv3[d] * 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][2] * func.vol;
        }

        std::vector<double> backAdv2L, backAdv3L;
        VecTool::resize(backAdv2L, main.dim);
        VecTool::resize(backAdv3L, main.dim);

        for(int d1=0; d1<main.dim; d1++){
            for(int d2=0; d2<main.dim; d2++){
                backAdv2L[d1] += adjoint.advk2[d2] * adjoint.dvk1dx[d1][d2];
                backAdv3L[d1] += adjoint.advk3[d2] * adjoint.dvk2dx[d1][d2];
            }
        }

        std::vector<double> backAdv2R, backAdv3R;
        VecTool::resize(backAdv2R, main.dim);
        VecTool::resize(backAdv3R, main.dim);

        for(int d1=0; d1<main.dim; d1++){
            for(int d2=0; d2<main.dim; d2++){
                backAdv2R[d1] += adjoint.advk2[d2] * adjoint.dvkdx[d1][d2];
                backAdv3R[d1] += adjoint.advk3[d2] * adjoint.dvk1dx[d1][d2];
            }
        }

        for(int d=0; d<main.dim; d++){
            value[p][0] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dwk2dx[d][0] * (backAdv3L[0] + backAdv3R[0]) * func.vol;
            value[p][0] += adjoint.tau * 1.5 * func.N[p] * adjoint.dwk1dx[d][0] * (backAdv2L[0] + backAdv2R[0]) * func.vol;
        }

        for(int d=0; d<main.dim; d++){
            value[p][1] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dwk2dx[d][1] * (backAdv3L[1] + backAdv3R[1]) * func.vol;
            value[p][1] += adjoint.tau * 1.5 * func.N[p] * adjoint.dwk1dx[d][1] * (backAdv2L[1] + backAdv2R[1]) * func.vol;
        }

        for(int d=0; d<main.dim; d++){
            value[p][2] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dwk2dx[d][2] * (backAdv3L[2] + backAdv3R[2]) * func.vol;
            value[p][2] += adjoint.tau * 1.5 * func.N[p] * adjoint.dwk1dx[d][2] * (backAdv2L[2] + backAdv2R[2]) * func.vol;
        }

        // SUPG pressure term
        for(int d=0; d<3; d++){
            value[p][0] += adjoint.tau * 1.5 * func.N[p] * adjoint.dpk1dx[d] * adjoint.dwk1dx[d][0] * func.vol;
            value[p][0] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dpk2dx[d] * adjoint.dwk2dx[d][0] * func.vol;
        }
        for(int d=0; d<3; d++){
            value[p][1] += adjoint.tau * 1.5 * func.N[p] * adjoint.dpk1dx[d] * adjoint.dwk1dx[d][1] * func.vol;
            value[p][1] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dpk2dx[d] * adjoint.dwk2dx[d][1] * func.vol;
        }
        for(int d=0; d<3; d++){
            value[p][2] += adjoint.tau * 1.5 * func.N[p] * adjoint.dpk1dx[d] * adjoint.dwk1dx[d][2] * func.vol;
            value[p][2] -= adjoint.tau * 0.5 * func.N[p] * adjoint.dpk2dx[d] * adjoint.dwk2dx[d][2] * func.vol;
        }

        // PSPG mass term 
        value[p][0] -= adjoint.tau * func.N[p] * adjoint.dqk1dx[0] / main.dt * func.vol;
        value[p][1] -= adjoint.tau * func.N[p] * adjoint.dqk1dx[1] / main.dt * func.vol;
        value[p][2] -= adjoint.tau * func.N[p] * adjoint.dqk1dx[2] / main.dt * func.vol;

        // PSPG advection term
        value[p][0] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][0] * adjoint.dqk1dx[0] * func.vol;
        value[p][0] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][0] * adjoint.dqk1dx[0] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][0] += adjoint.tau * 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * adjoint.dqk1dx[0] * func.vol;
        }
        value[p][0] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][0] * adjoint.dqk1dx[1] * func.vol;
        value[p][0] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][0] * adjoint.dqk1dx[1] * func.vol;
        value[p][0] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][0] * adjoint.dqk1dx[2] * func.vol;
        value[p][0] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][0] * adjoint.dqk1dx[2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][0] += adjoint.tau * 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][0] * adjoint.dqk2dx[d] * func.vol;
            value[p][0] += adjoint.tau * 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][0] * adjoint.dqk2dx[d] * func.vol;
        }

        value[p][1] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][1] * adjoint.dqk1dx[1] * func.vol;
        value[p][1] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][1] * adjoint.dqk1dx[1] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][1] += adjoint.tau * 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * adjoint.dqk1dx[1] * func.vol;
        }
        value[p][1] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][1] * adjoint.dqk1dx[0] * func.vol;
        value[p][1] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][1] * adjoint.dqk1dx[0] * func.vol;
        value[p][1] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][1] * adjoint.dqk1dx[2] * func.vol;
        value[p][1] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][1] * adjoint.dqk1dx[2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][1] += adjoint.tau * 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][1] * adjoint.dqk2dx[d] * func.vol;
            value[p][1] += adjoint.tau * 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][1] * adjoint.dqk2dx[d] * func.vol;
        }

        value[p][2] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[2][2] * adjoint.dqk1dx[2] * func.vol;
        value[p][2] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[2][2] * adjoint.dqk1dx[2] * func.vol;
        for(int d=0; d<main.dim; d++){
            value[p][2] += adjoint.tau * 0.5 * func.dNdx[p][d] * adjoint.advk2[d] * adjoint.dqk1dx[2] * func.vol;
        }
        value[p][2] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[0][2] * adjoint.dqk1dx[0] * func.vol;
        value[p][2] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[0][2] * adjoint.dqk1dx[0] * func.vol;
        value[p][2] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvk1dx[1][2] * adjoint.dqk1dx[1] * func.vol;
        value[p][2] += adjoint.tau * 0.5 * 1.5 * func.N[p] * adjoint.dvkdx[1][2] * adjoint.dqk1dx[1] * func.vol;   
        for(int d=0; d<main.dim; d++){
            value[p][2] += adjoint.tau * 0.5 * (-0.5) * func.N[p] * adjoint.dvk2dx[d][2] * adjoint.dqk2dx[d] * func.vol;
            value[p][2] += adjoint.tau * 0.5 * (-0.5) * func.N[p] * adjoint.dvk1dx[d][2] * adjoint.dqk2dx[d] * func.vol;
        }
    }
}

/***************************************************************
 * @brief Decide step length for both X and X0 at the same time.
 */
double InverseProblem::armijoCriteria(const double fk)
{
    const double c1 = 1e-2;
    double alpha = 1e0;

    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();

    std::vector<std::map<int, std::vector<double>>> vDirichletTmp;
    std::vector<std::map<int, std::vector<double>>> vDirichletNewTmp;
    vDirichletTmp.resize(main.timeMax);
    vDirichletNewTmp.resize(main.timeMax); 
    
    std::vector<std::vector<double>> v0Tmp;
    VecTool::resize(v0Tmp, main.grid.node.nNodesGlobal, main.dim);
    
    while(true){
        for(int in=0; in<main.grid.node.nNodesGlobal; in++){
            for(int d=0; d<main.dim; d++){
                v0Tmp[in][d] = X0[in][d] + alpha * (-gradInitVel[in][d]);
            }
        }
        
        for(int t=0; t<main.timeMax; t++){
            for(int ib=0; ib<n; ib++){
                int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
                std::vector<double> vecTmp(dim, 0e0);
                for(int d=0; d<dim; d++){
                    double value = X[t][ib][d] + alpha * (-grad[t][ib][d]);
                    vecTmp[d] = value;
                }
                vDirichletTmp[t][in] = vecTmp;
            }
        }
        for(int t=0; t<main.timeMax; t++){
            for(auto &pair : vDirichletTmp[t]){
                std::vector<double> vecTmp;
                int in = main.grid.node.mapNew[pair.first];
                for(auto &value : pair.second){
                    vecTmp.push_back(value);
                }
                vDirichletNewTmp[t][in] = vecTmp;
            }
        }

        main.solveUSNS(vDirichletNewTmp,
                       main.grid.dirichlet.pDirichletNew, v0Tmp);
        compCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int in=0; in<main.grid.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                tmp += -(gradInitVel[in][d] * gradInitVel[in][d]);
            }
        }
        for(int t=0; t<adjoint.timeMax; t++){
            if(t % main.snap.snapInterval == 0){
                for(int ib=0; ib<n; ib++){
                    for(int d=0; d<dim; d++){
                        tmp += -(grad[t][ib][d] * grad[t][ib][d]);
                    }
                }
            }
        }

        double l_tmp = fk + c1 * tmp * alpha;

        if(lk <= l_tmp){
            main.grid.node.v0 = v0Tmp;
            main.grid.dirichlet.vDirichlet = vDirichletTmp;
            main.grid.dirichlet.vDirichletNew = vDirichletNewTmp;
            break;
        }else{
            alpha = alpha * 5e-1;
            PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alpha);
        }
    }

    return alpha;
}


/***********************************
 * @brief Decide step length for X0.
 */
void InverseProblem::armijoCriteriaX0(const double fk)
{
    const double c1 = 1e-4;
    
    std::vector<std::vector<double>> v0Tmp;
    VecTool::resize(v0Tmp, main.grid.node.nNodesGlobal, main.dim);
    
    while(true){
        for(int in=0; in<main.grid.node.nNodesGlobal; in++){
            for(int d=0; d<main.dim; d++){
                v0Tmp[in][d] = X0[in][d] + alphaX0 * (-gradInitVel[in][d]);
            }
        }

        main.solveUSNS(main.grid.dirichlet.vDirichletNew, 
                       main.grid.dirichlet.pDirichletNew, v0Tmp);
        compCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int in=0; in<main.grid.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                tmp += -(gradInitVel[in][d] * gradInitVel[in][d]);
            }
        }

        double l_tmp = fk + c1 * tmp * alphaX0;

        if(lk <= l_tmp){
            main.grid.node.v0 = v0Tmp;
            break;
        }else{
            alphaX0 = alphaX0 * 5e-1;
            PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alphaX0);
        }
    }
}

/**********************************
 * @brief Decide step length for X.
 */
void InverseProblem::armijoCriteriaX(const double fk)
{
    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    const double c1 = 1e-4;

    std::vector<std::map<int, std::vector<double>>> vDirichletTmp;
    std::vector<std::map<int, std::vector<double>>> vDirichletNewTmp;
    vDirichletTmp.resize(main.timeMax);
    vDirichletNewTmp.resize(main.timeMax); 
 
    while(true){
        for(int t=0; t<main.timeMax; t++){
            for(int ib=0; ib<n; ib++){
                int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
                std::vector<double> vecTmp(dim, 0e0);
                for(int d=0; d<dim; d++){
                    double value = X[t][ib][d] + alphaX * (-grad[t][ib][d]);
                    vecTmp[d] = value;
                }
                vDirichletTmp[t][in] = vecTmp;
            }
        }
        for(int t=0; t<main.timeMax; t++){
            for(auto &pair : vDirichletTmp[t]){
                std::vector<double> vecTmp;
                int in = main.grid.node.mapNew[pair.first];
                for(auto &value : pair.second){
                    vecTmp.push_back(value);
                }
                vDirichletNewTmp[t][in] = vecTmp;
            }
        }

        main.solveUSNS(vDirichletNewTmp, main.grid.dirichlet.pDirichletNew, main.grid.node.v0);
        compCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int t=0; t<adjoint.timeMax; t++){
            if(t % main.snap.snapInterval == 0){
                for(int ib=0; ib<n; ib++){
                    for(int d=0; d<dim; d++){
                        tmp += -(grad[t][ib][d] * grad[t][ib][d]);
                    }
                }
            }
        }
        double l_tmp = fk + c1 * tmp * alphaX;

        if(lk <= l_tmp){
            main.grid.dirichlet.vDirichlet = vDirichletTmp;
            main.grid.dirichlet.vDirichletNew = vDirichletNewTmp;
            break;
        }else{
            alphaX = alphaX * 5e-1;
            PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alphaX);
        }
    }

}

/***********************************
 * @brief Decide step length for X0.
 */
double InverseProblem::armijoCriteriaX0_tmp(const double fk)
{
    if(isConverged_X0) return 0e0;

    const double c1 = 1e-3;
    double alpha = 1e0;
    
    std::vector<std::vector<double>> v0Tmp;
    VecTool::resize(v0Tmp, main.grid.node.nNodesGlobal, main.dim);
    
    int iterationCount = 0;
    const int maxIteration = 10;
    
    while(true){
        iterationCount++;
        if(iterationCount > maxIteration){
            PetscPrintf(MPI_COMM_WORLD, "Armijo X0 reached max iteration.\n");
            isConverged_X0 = true;
            break;
        }
        
        for(int in=0; in<main.grid.node.nNodesGlobal; in++){
            for(int d=0; d<main.dim; d++){
                v0Tmp[in][d] = X0[in][d] + alpha * (-gradInitVel[in][d]);
            }
        }

        main.solveUSNS(main.grid.dirichlet.vDirichletNew, 
                       main.grid.dirichlet.pDirichletNew, v0Tmp);
        compCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int in=0; in<main.grid.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                tmp += -(gradInitVel[in][d] * gradInitVel[in][d]);
            }
        }

        double l_tmp = fk + c1 * tmp * alpha;

        if(lk <= l_tmp){
            main.grid.node.v0 = v0Tmp;
            break;
        }else{
            alpha = alpha * 5e-1;
            PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alpha);
        }
    }

    return alpha;
}

/**********************************
 * @brief Decide step length for X.
 */
double InverseProblem::armijoCriteriaX_tmp(const double fk)
{
    if(isConverged_X) return 0e0;
    
    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    const double c1 = 1e-3;
    double alpha = 1e0;

    std::vector<std::map<int, std::vector<double>>> vDirichletTmp;
    std::vector<std::map<int, std::vector<double>>> vDirichletNewTmp;
    vDirichletTmp.resize(main.timeMax);
    vDirichletNewTmp.resize(main.timeMax); 

    int iterationCount = 0;
    const int maxIteration = 10;
 
    while(true){
        iterationCount++;
        if(iterationCount > maxIteration){
            PetscPrintf(MPI_COMM_WORLD, "Armijo X reached max iteration.\n");
            isConverged_X = true;
            break;
        }

        for(int t=0; t<main.timeMax; t++){
            for(int ib=0; ib<n; ib++){
                int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
                std::vector<double> vecTmp(dim, 0e0);
                for(int d=0; d<dim; d++){
                    double value = X[t][ib][d] + alpha * (-grad[t][ib][d]);
                    vecTmp[d] = value;
                }
                vDirichletTmp[t][in] = vecTmp;
            }
        }
        for(int t=0; t<main.timeMax; t++){
            for(auto &pair : vDirichletTmp[t]){
                std::vector<double> vecTmp;
                int in = main.grid.node.mapNew[pair.first];
                for(auto &value : pair.second){
                    vecTmp.push_back(value);
                }
                vDirichletNewTmp[t][in] = vecTmp;
            }
        }

        main.solveUSNS(vDirichletNewTmp, main.grid.dirichlet.pDirichletNew, main.grid.node.v0);
        compCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int t=0; t<adjoint.timeMax; t++){
            if(t % main.snap.snapInterval == 0){
                for(int ib=0; ib<n; ib++){
                    for(int d=0; d<dim; d++){
                        tmp += -(grad[t][ib][d] * grad[t][ib][d]);
                    }
                }
            }
        }
        double l_tmp = fk + c1 * tmp * alpha;

        if(lk <= l_tmp){
            main.grid.dirichlet.vDirichlet = vDirichletTmp;
            main.grid.dirichlet.vDirichletNew = vDirichletNewTmp;
            break;
        }else{
            alpha = alpha * 5e-1;
            PetscPrintf(MPI_COMM_WORLD, "Almijo %e %e %e %e\n", fk, lk, l_tmp, alpha);
        }
    }

    return alpha;
}