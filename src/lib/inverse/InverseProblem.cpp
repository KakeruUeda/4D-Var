#include "InverseProblem.h"

InverseProblem::InverseProblem(Config &conf):
app(conf.app),
main(conf), adjoint(conf), data(conf),
dim(conf.dim), outputDir(conf.outputDir), 
nOMP(conf.nOMP), timeMax(conf.timeMax),
rho(conf.rho), mu(conf.mu), dt(conf.dt),
alpha(conf.alpha), resistance(conf.resistance),
aCF(conf.aCF), bCF1(conf.bCF1), bCF2(conf.bCF2), gCF(conf.gCF),
loopMax(conf.loopMax), nControlNodesInCell(conf.nControlNodesInCell)
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
}

void InverseProblem::runSimulation()
{
    main.outputDomain();

    for(int loop=0; loop<loopMax; loop++){
        if(loop != 0) main.solveUSNS(app);
        calcCostFunction();
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
            vtiFile = outputDir + "/velocity/data" + to_string(t) + ".vti";
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
        }
    }

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

