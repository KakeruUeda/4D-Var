/**
 * @file InverseProblem.cpp
 * @author k.ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

InverseProblem::InverseProblem(Config &conf):
main(conf), adjoint(conf), data(conf), app(conf.app), vvox(conf.vvox), 
dim(conf.dim), outputDir(conf.outputDir), nOMP(conf.nOMP),
aCF(conf.aCF), bCF1(conf.bCF1), bCF2(conf.bCF2), gCF(conf.gCF),
loopMax(conf.loopMax), nControlNodesInCell(conf.nControlNodesInCell),
planeDir(conf.planeDir)
{
    std::string dir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = outputDir + "/adjoint";
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/feedback";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/solution";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/data";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}

void InverseProblem::runSimulation()
{
    main.outputDomain();
    guessInitialCondition();
    for(int loop=0; loop<loopMax; loop++){
        PetscPrintf(MPI_COMM_WORLD, "\nOpt itr. : %d\n", loop);
        compCostFunction();

        PetscPrintf(MPI_COMM_WORLD, "\ncostFunction = %e\n", costFunction.total);
        compFeedbackForce();
        compTimeInterpolatedFeedbackForce();
        output(loop);

        adjoint.solveAdjoint(main, outputDir, feedbackForceT, data.nData, loop);
        compOptimalCondition();
        double alpha = armijoCriteria(costFunction.total);
        PetscPrintf(MPI_COMM_WORLD, "\nalpha = %f\n", alpha);
        updataControlVariables(main, alpha);
    }
}

void InverseProblem::guessInitialCondition()
{
    main.compInitialCondition(main.grid.dirichlet.vDirichletNew, 
                              main.grid.dirichlet.pDirichletNew);
    main.solveUSNS(main.grid.dirichlet.vDirichletNew, 
                   main.grid.dirichlet.pDirichletNew);

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

    if(mpi.myId == 0){
        std::string vtuFile;
        vtuFile = main.outputDir + "/solution/velocityInitial.vtu";
        main.grid.output.exportSolutionVTU(vtuFile, main.grid.node, main.grid.cell, DataType::VELOCITY);  
    }
}

void InverseProblem::updataControlVariables(DirectProblem &main, const double alpha)
{   
    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    for(int t=0; t<main.timeMax; t++){
        for(int ib=0; ib<n; ib++){
            int in = adjoint.grid.dirichlet.controlBoundaryMap[ib];
            std::vector<double> vecTmp(dim, 0e0);
            for(int d=0; d<dim; d++){
                X[t][ib][d] += alpha * (-grad[t][ib][d]);
            }
        }
    }
}

void InverseProblem::output(const int loop)
{
    if(mpi.myId > 0) return;

    std::string vtiFile;
    std::string vtuFile;

    for(int t=0; t<main.snap.nSnapShot; t++){
        vtiFile = main.outputDir + "/data/data_" + to_string(loop) + "_" + to_string(t) + ".vti";
        main.grid.output.exportVelocityDataVTI(vtiFile, data, t, main.dim);
    }
    for(int t=0; t<main.timeMax; t++){
        vtuFile = outputDir + "/feedback/feedbackForceT_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        main.grid.output.exportFeedbackForceVTU(vtuFile, main.grid.node, main.grid.cell, feedbackForceT[t]);
        vtuFile = main.outputDir + "/solution/velocity_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        main.grid.output.exportMainVariablesVTU(vtuFile, main.grid.node, main.grid.cell, t, DataType::VELOCITY);
        vtuFile = main.outputDir + "/solution/pressure_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        main.grid.output.exportMainVariablesVTU(vtuFile, main.grid.node, main.grid.cell, t, DataType::PRESSURE);
    }

}

void InverseProblem::compCostFunction()
{  
    costFunction.term1 = 0e0;
    costFunction.term2 = 0e0;
    costFunction.term3 = 0e0;

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
                    data(ic).average(main.grid.cell, main.snap.v[t],
                             t, main.grid.cell.nNodesInCell, main.dim);
                }
            }
        break;

        case VoxelVelocity::INTERPOLATION:
            for(int t=0; t<main.snap.nSnapShot; t++){
                for(int ic=0; ic<data.nCellsGlobal; ic++){
                    data(ic).interpolate(main.grid.node, main.grid.cell, 
                                         main.snap.v[t], t, main.dim);
                }
            }
        break;
    
        default:
            PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
            exit(1);
    }

    // term1
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<data.nCellsGlobal; ic++){
            double dev = 0e0;
            double volume = data.dx * data.dy * data.dz;
            double deltaT = main.dt * main.snap.snapInterval;
            for(int d=0; d<dim; d++){
                data(ic).ve[t][d] = data(ic).vMRI[t][d] - data(ic).vCFD[t][d];
                dev += data(ic).ve[t][d] * data(ic).ve[t][d];
            }
            costFunction.term1 += 5e-1 * aCF * dev * volume * deltaT;
        }
    }

    int nc = adjoint.grid.dirichlet.nControlNodesInCell;

    Gauss gauss(2);
    Function func2d(nc, dim-1);

    // term2
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            for(int p=0; p<nc; p++){
                for(int d=0; d<dim-1; d++){
                    int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    func2d.xCurrent[p][d] = main.grid.node.x[n][planeDir[d]];
                }
            }
            double value = 0e0;
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    func2d.weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(func2d.N, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(func2d.dNdr, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralRegTerm2(func2d, value, ic, t);
                }
            }
            costFunction.term2 += 5e-1 * bCF1 * value;
        }
    }

    // term3
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            for(int p=0; p<nc; p++){
                for(int d=0; d<dim-1; d++){
                    int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    func2d.xCurrent[p][d] = main.grid.node.x[index][planeDir[d]];
                }
            }
            double value = 0e0;
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    func2d.weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(func2d.N, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(func2d.dNdr, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralRegTerm3(func2d, value, ic, t);
                }
            }
            costFunction.term3 += 5e-1 * bCF2 * value;
        }
    }
    costFunction.sum();
}

void InverseProblem::GaussIntegralRegTerm2(Function &func, double &value, const int ic, const int t)
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    double dxdr[2][2]; 

    MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
    func.detJ = MathCommon::compDeterminant_2x2(dxdr);
    func.vol = func.detJ * func.weight;

    std::vector<double> u(dim, 0e0);

    for(int p=0; p<nControlNodesInCell; p++){
        int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<dim; d++){
            u[d] += func.N[p] * main.snap.v[t][index][d];
        }
    }

    for(int d=0; d<dim; d++)
        value += u[d] * u[d] * func.vol;
}

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

    for(int p=0; p<nControlNodesInCell; p++){
        int n = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                dudx[d1][d2] += func.dNdx[p][d2] * main.snap.v[t][n][d1];
            }
        }
    }

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            value += dudx[d1][d2] * func.vol;
        }
    }
}

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

void InverseProblem::compTimeInterpolatedFeedbackForce()
{
    for(int in=0; in<adjoint.grid.node.nNodesGlobal; in++){
        std::vector<double> x, y1, y2, y3;
        for(int t=0; t<main.snap.nSnapShot; t++){
            double n = t * main.dt * main.snap.snapInterval;
            x.push_back(n);
            y1.push_back(feedbackForce[t][in][0]);
            y2.push_back(feedbackForce[t][in][1]);  
            y3.push_back(feedbackForce[t][in][2]); 
        }
        Spline splineX;
        Spline splineY;
        Spline splineZ;
        splineX.setPoints(x, y1);
        splineY.setPoints(x, y2);
        splineZ.setPoints(x, y3);
        for(int t=0; t<adjoint.timeMax; t++){
            double n = t * main.dt;
            if(t == 0){
                feedbackForceT[t][in][0] = feedbackForce[0][in][0];
                feedbackForceT[t][in][1] = feedbackForce[0][in][1];
                feedbackForceT[t][in][2] = feedbackForce[0][in][2];
            }else if(t == adjoint.timeMax - 1){
                feedbackForceT[t][in][0] = feedbackForce[main.snap.nSnapShot - 1][in][0];
                feedbackForceT[t][in][1] = feedbackForce[main.snap.nSnapShot - 1][in][1];
                feedbackForceT[t][in][2] = feedbackForce[main.snap.nSnapShot - 1][in][2];
            }else{
                feedbackForceT[t][in][0] = splineX(n);
                feedbackForceT[t][in][1] = splineY(n);
                feedbackForceT[t][in][2] = splineZ(n);
            }
        }
    }
}


void InverseProblem::compOptimalCondition()
{
    int nc = adjoint.grid.dirichlet.nControlNodesInCell;
    Function func2d(nc, dim-1);
    Gauss gauss(2);

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

        // term1
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            double value[4][3];
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
                    GaussIntegralOptimalConditionTerm1(func2d, value, ic, t);
                }
            }
            for(int p=0; p<nc; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += bCF1 * value[p][d];
                }
            }
        }

        // term2
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            double value[4][3];
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
                    GaussIntegralOptimalConditionTerm2(func2d, value, ic, t);
                }
            }
            for(int p=0; p<nc; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += bCF2 * value[p][d];
                }
            }
        }

        // term3
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            double value[4][3];
            for(int p=0; p<nc; p++){
                for(int d=0; d<3; d++){
                    value[p][d] = 0e0;
                }
            }
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
                    GaussIntegralOptimalConditionTerm3(func2d, value, ic, t);
                }
            }
            for(int p=0; p<nControlNodesInCell; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += value[p][d];
                }
            }
        }

        for(int in=0; in<adjoint.grid.node.nNodesGlobal; in++){
            if(adjoint.grid.dirichlet.isBoundaryEdge[in]){
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] = 0e0;
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

    if(mpi.myId == 0){
        for(int t=0; t<main.timeMax; t++){
            std::string vtuFile;
            vtuFile = outputDir + "/feedback/gradWholeNode_" + to_string(t) + ".vtu";
            main.grid.output.exportFeedbackForceVTU(vtuFile, main.grid.node, main.grid.cell, gradWholeNode[t]);
        }
    }

}

void InverseProblem::GaussIntegralOptimalConditionTerm1(Function &func, double (&value)[4][3], const int ic, const int t)
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

void InverseProblem::GaussIntegralOptimalConditionTerm2(Function &func, double (&value)[4][3], const int ic, const int t)
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

void InverseProblem::GaussIntegralOptimalConditionTerm3(Function &func, double (&value)[4][3], const int ic, const int t)
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

double InverseProblem::armijoCriteria(const double fk)
{
    int n = adjoint.grid.dirichlet.controlBoundaryMap.size();
    const double c1 = 1e-2;
    double alpha = 1e0;

    std::vector<std::map<int, std::vector<double>>> vDirichletTmp;
    std::vector<std::map<int, std::vector<double>>> vDirichletNewTmp;
    vDirichletTmp.resize(adjoint.timeMax);
    vDirichletNewTmp.resize(adjoint.timeMax); 
    
    std::vector<std::map<int, double>> pDirichletNewTmp;
    pDirichletNewTmp = adjoint.grid.dirichlet.pDirichletNew;
 
    while(10){
        for(int t=0; t<adjoint.timeMax; t++){
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
        for(int t=0; t<adjoint.timeMax; t++){
            for(auto &pair : vDirichletTmp[t]){
                std::vector<double> vecTmp;
                int in = adjoint.grid.node.mapNew[pair.first];
                for(auto &value : pair.second){
                    vecTmp.push_back(value);
                }
                vDirichletNewTmp[t][in] = vecTmp;
            }
        }

        main.compInitialCondition(vDirichletNewTmp, pDirichletNewTmp);
        main.solveUSNS(vDirichletNewTmp, pDirichletNewTmp);
        compCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int t=0; t<adjoint.timeMax; t++){
            //if(t % main.snap.snapInterval == 0){
                for(int ib=0; ib<n; ib++){
                    for(int d=0; d<dim; d++){
                        tmp += -(grad[t][ib][d] * grad[t][ib][d]);
                    }
                }
            //}
        }
        double l_tmp = fk + c1 * tmp * alpha;

        // printf("Almijo %e %e %e\n",fk,lk,l_tmp);

        if(lk <= l_tmp){
            main.grid.dirichlet.vDirichlet = vDirichletTmp;
            main.grid.dirichlet.vDirichletNew = vDirichletNewTmp;
            break;
        }else{
            alpha = alpha * 5e-1;
            printf("Almijo %e %e %e %e\n", fk, lk, l_tmp, alpha);
        }
    }

    return alpha;
}