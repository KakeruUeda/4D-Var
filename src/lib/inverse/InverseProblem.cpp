#include "InverseProblem.h"

InverseProblem::InverseProblem(Config &conf):
main(conf), adjoint(conf), data(conf), app(conf.app), 
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

void InverseProblem::initialize(Config &conf)
{
    PetscPrintf(MPI_COMM_WORLD, "\n*** Main initialize ***\n\n");

    main.nu = main.mu / main.rho;
    main.Re = 1e0 / main.nu; 

    VecTool::resize(main.grid.node.v, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.vPrev, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.vt, main.timeMax, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.p, main.grid.node.nNodesGlobal);
    VecTool::resize(main.grid.node.pt, main.timeMax, main.grid.node.nNodesGlobal);
    VecTool::resize(main.snap.v, main.snap.nSnapShot, main.grid.nNodesGlobal, dim);

    main.grid.dirichlet.initialize(conf);
    main.grid.cell.initialize(conf);
    main.grid.node.initialize(conf);
    main.grid.prepareMatrix(main.petsc, main.outputDir, main.timeMax);

    VecTool::resize(main.petsc.solution, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValue, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueNew, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueInit, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueNewInit, main.grid.nDofsGlobal);

    PetscPrintf(MPI_COMM_WORLD, "\n*** Adjoint initialize ***\n\n");

    adjoint.nu = adjoint.mu / adjoint.rho;
    adjoint.Re = 1e0 / adjoint.nu; 

    VecTool::resize(adjoint.grid.node.w, adjoint.grid.node.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.q, adjoint.grid.node.nNodesGlobal);
    VecTool::resize(adjoint.grid.node.l, main.grid.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.lt, adjoint.timeMax, adjoint.grid.nNodesGlobal, dim);

    adjoint.grid.dirichlet.initializeAdjoint(conf);
    adjoint.grid.cell.initializeAdjoint(conf);
    adjoint.grid.node.initializeAdjoint(conf, adjoint.grid.dirichlet.controlBoundaryMap);
    adjoint.grid.prepareMatrix(adjoint.petsc, outputDir, adjoint.timeMax);

    VecTool::resize(adjoint.petsc.solution, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValue, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValueNew, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValueInit, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValueNewInit, adjoint.grid.nDofsGlobal);
    VecTool::resize(feedbackForce, main.snap.nSnapShot, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(gradWholeNode, main.timeMax, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(grad, main.timeMax, adjoint.grid.dirichlet.controlBoundaryMap.size(), dim);
    VecTool::resize(X, main.timeMax, adjoint.grid.dirichlet.controlBoundaryMap.size(), dim);

    data.initialize(conf, main.grid.node, main.grid.cell, dim);
}

void InverseProblem::runSimulation()
{
    main.outputDomain();

    for(int loop=0; loop<loopMax; loop++){
        if(loop != 0) main.solveUSNS(app);
        calcCostFunction();
        calcFeedbackForce();
        output(loop);
        
        adjoint.solveAdjointEquation(main, outputDir, feedbackForce, data.nData, loop);
        calcOptimalCondition();

        double alpha;
        alpha = armijoCriteria(costFunction.total);
        PetscPrintf(MPI_COMM_WORLD, "\nalpha = %f\n", alpha);
    }
}

void InverseProblem::output(const int loop)
{
    if(mpi.myId > 0) return;

    for(int t=0; t<main.snap.nSnapShot; t++){
        std::string vtuFile;
        vtuFile = outputDir + "/feedback/feedback_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        main.grid.output.exportFeedbackForceVTU(vtuFile, main.grid.node, main.grid.cell, feedbackForce[t]);
    }
}

void InverseProblem::calcCostFunction()
{  
    costFunction.term1 = 0e0;
    costFunction.term2 = 0e0;
    costFunction.term3 = 0e0;

    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<data.nCellsGlobal; ic++){
            data(ic).averageVelocity(main.grid.cell, main.snap.v[t],
                           t, main.grid.cell.nNodesInCell, main.dim);
        }
    }

    for(int t=0; t<main.snap.nSnapShot; t++){
        if(mpi.myId == 0){
            std::string vtiFile;
            vtiFile = outputDir + "/data/data" + to_string(t) + ".vti";
            main.grid.output.exportVelocityDataVTI(vtiFile, data, t, main.dim);
        }
    }
    
    // term1
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<data.nCellsGlobal; ic++){
            double dev = 0e0;
            double volume = data.dx * data.dy * data.dz;
            for(int d=0; d<dim; d++){
                data(ic).ve[t][d] = data(ic).vCFD[t][d] - data(ic).vMRI[t][d];
                dev += data(ic).ve[t][d] * data(ic).ve[t][d];
            }
            costFunction.term1 += 5e-1 * aCF * dev * volume;
        }
    }

    std::vector<double> N2D;
    std::vector<std::vector<double>> dNdr2D;
    std::vector<std::vector<double>> xCurrent2D;

    N2D.resize(nControlNodesInCell, 0e0);
    dNdr2D.resize(nControlNodesInCell, std::vector<double>(dim-1, 0e0));
    xCurrent2D.resize(nControlNodesInCell, std::vector<double>(dim-1, 0e0));

    Gauss gauss(2);
    int planeDir[2] = {1, 2};
    double value;

    // term2
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            for(int p=0; p<nControlNodesInCell; p++){
                for(int d=0; d<dim-1; d++){
                    int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    xCurrent2D[p][d] = main.grid.node.x[index][planeDir[d]];
                }
            }
            double value = 0e0;
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    double weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralRegTerm1(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            costFunction.term2 += 5e-1 * bCF1 * value;
        }
    }

    // term3
    for(int t=0; t<main.snap.nSnapShot; t++){
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            for(int p=0; p<nControlNodesInCell; p++){
                for(int d=0; d<dim-1; d++){
                    int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    xCurrent2D[p][d] = main.grid.node.x[index][planeDir[d]];
                }
            }
            double value = 0e0;
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    double weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralRegTerm2(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            costFunction.term3 += 5e-1 * bCF2 * value;
        }
    }
    costFunction.sum();
}

void InverseProblem::GaussIntegralRegTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                           std::vector<std::vector<double>> &xCurrent, double &value, 
                                           const double weight, const int ic, const int t)
{
    double dxdr[2][2]; 
    MathFEM::calc_dxdr2D(dxdr, dNdr, xCurrent, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

    double u[3] = {0e0, 0e0, 0e0};

    for(int p=0; p<nControlNodesInCell; p++){
        int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<dim; d++){
            u[d] += N[p] * main.snap.v[t][index][d];
        }
    }

    for(int d=0; d<dim; d++)
        value += u[d] * u[d] * detJ * weight;
}

void InverseProblem::GaussIntegralRegTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                           std::vector<std::vector<double>> &xCurrent, double &value, 
                                           const double weight, const int ic, const int t)
{
    double dxdr[2][2]; 
    std::vector<std::vector<double>> dNdx;
    dNdx.resize(nControlNodesInCell, std::vector<double>(dim-1, 0e0));

    MathFEM::calc_dxdr2D(dxdr, dNdr, xCurrent, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
    MathFEM::calc_dNdx2D(dNdx, dNdr, dxdr, nControlNodesInCell);

    double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

    for(int p=0; p<nControlNodesInCell; p++){
        int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                dudx[d1][d2] = dNdx[p][d2] * main.snap.v[t][index][d1];
            }
        }
    }

    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<2; d2++){
            value += dudx[d1][d2] * detJ * weight;
        }
    }
}

void InverseProblem::calcFeedbackForce()
{   
    for(int t=0; t<main.snap.nSnapShot; t++){
        int n1 = main.grid.cell.nCellsGlobal;
        int n2 = main.grid.cell.nNodesInCell;

        std::vector<std::vector<std::vector<std::vector<double>>>> vEX;
        VecTool::resize(vEX, data.nz+2, data.ny+2, data.nx+2, dim);
        calcEdgeValue(vEX, t);

        for(int ic=0; ic<main.grid.cell.nCellsGlobal; ic++){
            std::vector<std::vector<double>> xCurrent;
            xCurrent.resize(main.grid.cell.nNodesInCell, std::vector<double>(dim, 0e0));

            for(int p=0; p<main.grid.cell.nNodesInCell; p++){
                for(int d=0; d<dim; d++){
                    xCurrent[p][d] = main.grid.cell(ic).x[p][d];
                }
            }

            std::vector<double> N;
            std::vector<std::vector<double>> dNdr;
            N.resize(main.grid.cell.nNodesInCell, 0e0);
            dNdr.resize(main.grid.cell.nNodesInCell, std::vector<double>(dim, 0e0));
            int nGaussPoint = 2;
            Gauss gauss(nGaussPoint);

            for(int i1=0; i1<nGaussPoint; i1++){
                for(int i2=0; i2<nGaussPoint; i2++){
                    for(int i3=0; i3<nGaussPoint; i3++){
                        double dxdr[3][3];
                        ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                        ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                        MathFEM::calc_dxdr(dxdr, dNdr, xCurrent, main.grid.cell.nNodesInCell);
                        double detJ = MathCommon::calcDeterminant_3x3(dxdr);
                        double weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                        double feedback[3] = {0e0, 0e0, 0e0};
                        double point[3] = {0e0, 0e0, 0e0};
                        for(int d=0; d<dim; d++){
                            for(int p=0; p<main.grid.cell.nNodesInCell; p++) 
                                point[d] += N[p] * xCurrent[p][d];
                        }
                        calcInterpolatedFeeback(xCurrent, feedback, vEX, point);
                        feedbackGaussIntegral(N, feedback, detJ, weight, ic, t);
                    }
                }
            }
        }
    }
}

void InverseProblem::calcEdgeValue(std::vector<std::vector<std::vector<std::vector<double>>>> &vEX, const int t)
{
    for(int k=0; k<data.nz+2; k++){
        for(int j=0; j<data.ny+2; j++){
            for(int i=0; i<data.nx+2; i++){
                for(int d=0; d<dim; d++){
                    vEX[k][j][i][d] = 0e0;
                }
            }
        }
    }
    for(int k=0; k<data.nz; k++){
        for(int j=0; j<data.ny; j++){
            for(int i=0; i<data.nx; i++){
                for(int d=0; d<dim; d++){
                    vEX[k+1][j+1][i+1][d] 
                    = data(k, j, i).vMRI[t][d];
                }
            }
        }
    }

}

void InverseProblem::calcInterpolatedFeeback(std::vector<std::vector<double>> &xCurrent, double (&feedback)[3], 
                                             std::vector<std::vector<std::vector<std::vector<double>>>> &vEX, 
                                             double (&point)[3])
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

    std::vector<double> N(main.grid.cell.nNodesInCell, 0e0);
    ShapeFunction3D::C3D8_N(N, s, t, u);

    for(int d=0; d<dim; d++){
        feedback[d] = N[0]*vEX[iz][iy][ix][d]       +  N[1]*vEX[iz][iy][ix+1][d]
                    + N[2]*vEX[iz][iy+1][ix+1][d]   +  N[3]*vEX[iz][iy+1][ix][d]
                    + N[4]*vEX[iz+1][iy][ix][d]     +  N[5]*vEX[iz+1][iy][ix+1][d]
                    + N[6]*vEX[iz+1][iy+1][ix+1][d] +  N[7]*vEX[iz+1][iy+1][ix][d];
    }
}

void InverseProblem::feedbackGaussIntegral(std::vector<double> &N, double (&feedback)[3], 
                                           const double detJ, const double weight, const int ic, const int t)
{
    for(int d=0; d<dim; d++){
        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            int in = main.grid.cell(ic).node[p];
            feedbackForce[t][in][d] += aCF * N[p] * feedback[d] * detJ * weight;
        }
    }
}


void InverseProblem::calcFeedbackForce2()
{
    for(int t=0; t<main.snap.nSnapShot; t++){
        int n1 = main.grid.cell.nCellsGlobal;
        int n2 = main.grid.cell.nNodesInCell;

        std::vector<std::vector<double>> velCurrent;
        std::vector<std::vector<double>> xCurrent;
        velCurrent.resize(main.grid.cell.nNodesInCell, std::vector<double>(dim, 0e0));
        xCurrent.resize(main.grid.cell.nNodesInCell, std::vector<double>(dim, 0e0));

        for(int ic=0; ic<data.nCellsGlobal; ic++){
            for(int ic2=0; ic2<data(ic).cellChildren.size(); ic2++){
                int voxelId = ic;
                int cellId = data(ic).cellChildren[ic2];
                for(int p=0; p<main.grid.cell.nNodesInCell; p++){
                    for(int d=0; d<dim; d++){
                        velCurrent[p][d] = main.snap.v[t][main.grid.cell(cellId).node[p]][d];
                        xCurrent[p][d] = main.grid.cell(cellId).x[p][d];
                    }
                }
                std::vector<double> N;
                std::vector<std::vector<double>> dNdr;
                N.resize(main.grid.cell.nNodesInCell, 0e0);
                dNdr.resize(main.grid.cell.nNodesInCell, std::vector<double>(dim, 0e0));
                int nGaussPoint = 2;
                Gauss gauss(nGaussPoint);
                double detJ, weight;
                for(int i1=0; i1<nGaussPoint; i1++){
                    for(int i2=0; i2<nGaussPoint; i2++){
                        for(int i3=0; i3<nGaussPoint; i3++){
                            double dxdr[3][3];
                            ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                            ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                            MathFEM::calc_dxdr(dxdr, dNdr, xCurrent, main.grid.cell.nNodesInCell);
                            detJ = MathCommon::calcDeterminant_3x3(dxdr);
                            weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                            feedbackGaussIntegral2(N, xCurrent, velCurrent, detJ, weight, voxelId, cellId, t);
                        }
                    }
                }
            }
        }
    }
}

void InverseProblem::feedbackGaussIntegral2(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent,
                                           std::vector<std::vector<double>> &velCurrent, const double detJ,
                                           const double weight, const int voxelId, const int cellId, const int t)
{
    for(int d=0; d<dim; d++){
        for(int p=0; p<main.grid.cell.nNodesInCell; p++){
            int index = main.grid.cell(cellId).node[p];
            feedbackForce[t][index][d] += N[p] * (data(voxelId).vMRI[t][d] - velCurrent[p][d]) * detJ * weight;
        }
    }
}

void InverseProblem::calcOptimalCondition()
{
    std::vector<double> N2D;
    std::vector<std::vector<double>> dNdr2D;
    std::vector<std::vector<double>> xCurrent2D;

    N2D.resize(nControlNodesInCell, 0e0);
    dNdr2D.resize(nControlNodesInCell, std::vector<double>(dim-1, 0e0));
    xCurrent2D.resize(nControlNodesInCell, std::vector<double>(dim-1, 0e0));

    Gauss gauss(2);

    for(int t=0; t<main.timeMax; t++){
        // term1
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            double value[4][3];
            for(int p=0; p<nControlNodesInCell; p++){
                for(int d=0; d<dim-1; d++){
                    int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    xCurrent2D[p][d] = main.grid.node.x[in][planeDir[d]];
                }
            }
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    double weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralOptimalConditionTerm1(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            for(int p=0; p<nControlNodesInCell; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += bCF1 * value[p][d];
                }
            }
        }

        // term2
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            double value[4][3];
            for(int p=0; p<nControlNodesInCell; p++){
                for(int d=0; d<dim-1; d++){
                    int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    xCurrent2D[p][d] = main.grid.node.x[in][planeDir[d]];
                }
            }
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    double weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralOptimalConditionTerm2(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            for(int p=0; p<nControlNodesInCell; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += bCF2 * value[p][d];
                }
            }
        }

        // term3
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
            double value[4][3];
            for(int p=0; p<nControlNodesInCell; p++){
                for(int d=0; d<dim-1; d++){
                    int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                    xCurrent2D[p][d] = main.grid.node.x[in][planeDir[d]];
                }
            }
            for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                    double weight = gauss.weight[i1] * gauss.weight[i2];
                    ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
                    ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
                    GaussIntegralOptimalConditionTerm3(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            for(int p=0; p<nControlNodesInCell; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim; d++){
                    gradWholeNode[t][in][d] += value[p][d];
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

}

void InverseProblem::GaussIntegralOptimalConditionTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                         std::vector<std::vector<double>> &xCurrent, 
                                                         double (&value)[4][3], const double weight, 
                                                         const int ic, const int t)
{
    double vel[3] = {0e0, 0e0, 0e0};
    double dxdr[2][2]; 
    MathFEM::calc_dxdr2D(dxdr, dNdr, xCurrent, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

    for(int p=0; p<nControlNodesInCell; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<3; d++){
            vel[d] += N[p] * main.grid.node.vt[t][in][d];
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] = 0e0;
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] += vel[d] * N[p] * detJ * weight;
        }
    }
}

void InverseProblem::GaussIntegralOptimalConditionTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                         std::vector<std::vector<double>> &xCurrent, 
                                                         double (&value)[4][3], const double weight, 
                                                         const int ic, const int t)
{
    double dxdr[2][2]; 
    double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    std::vector<std::vector<double>> dNdx;
    dNdx.resize(nControlNodesInCell, std::vector<double>(2, 0e0));

    MathFEM::calc_dxdr2D(dxdr, dNdr, xCurrent, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
    MathFEM::calc_dNdx2D(dNdx, dNdr, dxdr, nControlNodesInCell);

    for(int p=0; p<nControlNodesInCell; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                dudx[d1][d2] = dNdx[p][d2] * main.grid.node.vt[t][in][d1];
            }
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] = 0e0;
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                value[p][d1] += dudx[d1][d2] * N[p]* detJ * weight;
            }
        }
    }
}

void InverseProblem::GaussIntegralOptimalConditionTerm3(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                         std::vector<std::vector<double>> &xCurrent, 
                                                         double (&value)[4][3], const double weight, 
                                                         const int ic, const int t)
{
    double gpValue[3] = {0e0, 0e0, 0e0};
    double dxdr[2][2]; 
    MathFEM::calc_dxdr2D(dxdr, dNdr, xCurrent, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

    for(int p=0; p<nControlNodesInCell; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<3; d++){
            gpValue[d] -= N[p] * adjoint.grid.node.lt[t][in][d];
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] = 0e0;
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] += gpValue[d] * N[p] * detJ * weight;
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

        main.solveUSNS(vDirichletNewTmp, pDirichletNewTmp);
        calcCostFunction();

        double lk = costFunction.total;

        double tmp = 0e0;
        for(int t=0; t<adjoint.timeMax; t++){
            for(int ib=0; ib<n; ib++){
                for(int d=0; d<dim; d++){
                    tmp += -(grad[t][ib][d] * grad[t][ib][d]);
                }
            }
        }
        double l_tmp = fk + c1 * tmp * alpha;

        // printf("Almijo %e %e %e\n",fk,lk,l_tmp);

        if (lk <= l_tmp){
            break;
        }else{
            alpha = alpha * 5e-1;
            printf("Almijo %e %e %e %e\n", fk, lk, l_tmp, alpha);
        }
    }

    return alpha;
}