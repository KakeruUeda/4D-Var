#include "InverseProblem.h"

InverseProblem::InverseProblem(Config &conf):
app(conf.app), main(conf), adjoint(conf), data(conf),
dim(conf.dim), outputDir(conf.outputDir), 
nOMP(conf.nOMP), timeMax(conf.timeMax),
rho(conf.rho), mu(conf.mu), dt(conf.dt),
alpha(conf.alpha), resistance(conf.resistance),
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
    dir = outputDir + "/velocity";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/pressure";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/lambda";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/domain";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}

void InverseProblem::initialize(Config &conf)
{
    PetscPrintf(MPI_COMM_WORLD, "\n*** Main initialize ***\n\n");
    main.grid.node.v.resize(main.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    main.grid.node.vPrev.resize(main.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    main.grid.node.vt.resize(timeMax, std::vector<std::vector<double>>(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    main.grid.node.p.resize(main.grid.node.nNodesGlobal, 0e0);
    main.grid.node.pt.resize(timeMax, std::vector<double>(main.grid.node.nNodesGlobal, 0e0));
    main.snap.v.resize(main.snap.nSnapShot, std::vector<std::vector<double>>(main.grid.nNodesGlobal, std::vector<double>(dim, 0e0)));
    main.grid.dirichlet.initialize(conf);
    main.grid.cell.initialize(conf);
    main.grid.node.initialize(conf);
    main.grid.prepareMatrix(main.petsc, main.outputDir);

    PetscPrintf(MPI_COMM_WORLD, "\n*** Adjoint initialize ***\n\n");
    adjoint.grid.node.v.resize(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    adjoint.grid.node.vPrev.resize(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    adjoint.grid.node.vt.resize(timeMax, std::vector<std::vector<double>>(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    adjoint.grid.node.p.resize(adjoint.grid.node.nNodesGlobal, 0e0);
    adjoint.grid.node.pt.resize(timeMax, std::vector<double>(adjoint.grid.node.nNodesGlobal, 0e0));
    adjoint.grid.node.lambda.resize(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    adjoint.grid.node.lambdat.resize(timeMax, std::vector<std::vector<double>>(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    adjoint.grid.dirichlet.initializeAdjoint(conf);
    adjoint.grid.cell.initializeAdjoint(conf);
    adjoint.grid.node.initializeAdjoint(conf, adjoint.grid.dirichlet.controlBoundaryMap);
    adjoint.grid.prepareMatrix(adjoint.petsc, outputDir);

    data.initialize(conf, main.grid.node, main.grid.cell, dim);

    feedbackForce.resize(main.snap.nSnapShot, std::vector<std::vector<double>>(main.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    gradWholeNode.resize(main.timeMax, std::vector<std::vector<double>>(main.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    grad.resize(main.timeMax, std::vector<std::vector<double>>(adjoint.grid.dirichlet.controlBoundaryMap.size(), std::vector<double>(dim, 0e0)));
}

void InverseProblem::runSimulation()
{
    main.outputDomain();

    for(int loop=0; loop<loopMax; loop++){
        if(loop != 0) main.solveUSNS(app);
        calcCostFunction();
        calcFeedbackForce();
        adjoint.solveAdjointEquation(main, outputDir, feedbackForce);
        calcOptimalCondition();
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
            vtiFile = outputDir + "/domain/data" + to_string(t) + ".vti";
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
            costFunction.term1 += 5e-1 * alpha * dev * volume;
        }
    }

    std::vector<double> N2D;
    std::vector<std::vector<double>> dNdr2D;
    std::vector<std::vector<double>> xCurrent2D;

    N2D.resize(nControlNodesInCell, 0e470);
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
    std::vector<std::vector<double>> dNdx;
    dNdx.resize(nControlNodesInCell, std::vector<double>(dim-1, 0e0));

    MathFEM::calc_dxdr2D(dxdr, dNdr, xCurrent, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
    MathFEM::calc_dNdx2D(dNdx, dNdr, dxdr, nControlNodesInCell);

    double u[dim] = {0e0, 0e0, 0e0};

    for(int p=0; p<nControlNodesInCell; p++){
        int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<dim; d++){
            u[d] = main.snap.v[t][index][d];
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

    double dudx[dim][dim-1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

    for(int p=0; p<nControlNodesInCell; p++){
        int index = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<dim; d1++){
            for(int d2=0; d2<dim-1; d2++){
                dudx[d1][d2] = dNdx[p][d2] * main.snap.v[t][index][d1];
            }
        }
    }

    for(int d1=0; d1<dim; d1++){
        for(int d2=0; d2<dim-1; d2++){
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
        vEX.resize(data.nz+2, std::vector<std::vector<std::vector<double>>>
                  (data.ny+2, std::vector<std::vector<double>>
                  (data.nx+2, std::vector<double>(dim, 0e0))));
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
    
    for(int t=0; t<main.snap.nSnapShot; t++){
        if(mpi.myId == 0){
            std::string vtuFile;
            vtuFile = outputDir + "/domain/force" + to_string(t) + ".vtu";
            main.grid.output.exportAdjointVTU(vtuFile, main.grid.node, main.grid.cell, feedbackForce[t], DataType::FEEDBACK);
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
            int np = main.grid.cell(ic).node[p];
            feedbackForce[t][np][d] += aCF * N[p] * feedback[d] * detJ * weight;
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
    
    for(int t=0; t<main.snap.nSnapShot; t++){
        if(mpi.myId == 0){
            std::string vtuFile;
            vtuFile = outputDir + "/domain/force" + to_string(t) + ".vtu";
            main.grid.output.exportAdjointVTU(vtuFile, main.grid.node, main.grid.cell, feedbackForce[t], DataType::FEEDBACK);
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
    double dxdr[2][2];
    double value[4][3];

    for(int t=0; t<main.timeMax; t++){
        // term1
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
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
                    GaussIntegralOpttimalConditionTerm1(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            for(int p=0; p<nControlNodesInCell; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim-1; d++){
                    gradWholeNode[t][in][d] += bCF1 * value[p][d];
                }
            }
        }

        // term1
        for(int ic=0; ic<adjoint.grid.dirichlet.controlNodeInCell.size(); ic++){
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
                    GaussIntegralOpttimalConditionTerm2(N2D, dNdr2D, xCurrent2D, value, weight, ic, t);
                }
            }
            for(int p=0; p<nControlNodesInCell; p++){
                int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
                for(int d=0; d<dim-1; d++){
                    gradWholeNode[t][in][d] += bCF1 * value[p][d];
                }
            }
        }
    }

}

void InverseProblem::GaussIntegralOpttimalConditionTerm1(std::vector<double> &N2D, std::vector<std::vector<double>> &dNdr2D, 
                                                         std::vector<std::vector<double>> &xCurrent2D, 
                                                         double (&value)[4][3], const double weight, 
                                                         const int ic, const int t)
{
    double vel[3] = {0e0, 0e0, 0e0};
    double dxdr[2][2]; 
    MathFEM::calc_dxdr2D(dxdr, dNdr2D, xCurrent2D, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

    for(int p=0; p<nControlNodesInCell; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d=0; d<3; d++){
            vel[d] += N2D[p] * main.grid.node.vt[t][in][planeDir[d]];
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] = 0e0;
        }
    }
    for(int p=0; p<nControlNodesInCell; p++){
        for(int d=0; d<3; d++){
            value[p][d] += vel[d] * N2D[p] * detJ * weight;
        }
    }
}

void InverseProblem::GaussIntegralOpttimalConditionTerm2(std::vector<double> &N2D, std::vector<std::vector<double>> &dNdr2D, 
                                                         std::vector<std::vector<double>> &xCurrent2D, 
                                                         double (&value)[4][3], const double weight, 
                                                         const int ic, const int t)
{
    double dxdr[2][2]; 
    double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    std::vector<std::vector<double>> dNdx2D;
    dNdx2D.resize(nControlNodesInCell, std::vector<double>(2, 0e0));

    MathFEM::calc_dxdr2D(dxdr, dNdr2D, xCurrent2D, nControlNodesInCell);
    double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
    MathFEM::calc_dNdx2D(dNdx2D, dNdr2D, dxdr, nControlNodesInCell);

    for(int p=0; p<nControlNodesInCell; p++){
        int in = adjoint.grid.dirichlet.controlNodeInCell[ic][p];
        for(int d1=0; d1<3; d1++){
            for(int d2=0; d2<2; d2++){
                dudx[d1][d2] = dNdx2D[p][d2] * main.grid.node.vt[t][in][planeDir[d1]];
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
                value[p][d1] += dudx[d1][d2] * N2D[p]* detJ * weight;
            }
        }
    }
}

