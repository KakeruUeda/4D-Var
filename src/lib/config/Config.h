/**
 * @file Config.cpp
 * @author K.Ueda
 * @date May, 2024
*/

#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <fstream>
#include <sstream>
#include "TextParser.h"
#include "Array.h"
#include "MyMPI.h"
#include "Define.h"
#include "Tool.h"

extern MyMPI mpi;

enum class Application
{
    STRGRID = 0, SNS = 1, USNS = 2,
    VOXELDATA = 3,  FDVAR = 4, 
    FLOWRATE = 5, MAE = 6
};

enum class GridType
{
    STRUCTURED = 0, 
    UNSTRUCTURED = 1
};

enum class ControlBoundary
{
    left = 0, right = 1, top = 2, bottom = 3,
    front = 4, back = 5
};

enum class VoxelVelocity
{
    POINTSPREAD = 0,
    AVERAGE = 1,
    INTERPOLATION = 2
};

enum class CrossSection
{
    XY = 0,
    YZ = 1,
    ZX = 2
};

class Config
{
    public:
        Config(std::string inputFile, std::string appName);
        ~Config(){};

        TextParser tp;
        Application app;
        GridType gridType;
        ControlBoundary controlBoundary;
        VoxelVelocity vvox;
        CrossSection crossSection;

        // Basic parameter
        int dim, nOMP;
        int outputItr;
        std::string outputDir;

        // Physical parameter
        double rho, mu, Re;

        double NRtolerance;

        // CostFunction parameter
        double aCF, bCF, gCF;
        int loopMax;

        // Time parameter
        double dt;
        int timeMax;
        int pulsatileFlow;
        int pulseBeginItr;
        double T;

        // Darcy Parameter
        double alpha, resistance;
    
        // Grid parameter
        int extractFluid;
        int nx, ny, nz;
        double lx, ly, lz;
        double dx, dy, dz;
        int nxNodes, nyNodes, nzNodes;
        int nxCells, nyCells, nzCells;
        int nCellsGlobal, nNodesGlobal, nNodesInCell;

        // SnapShot parameter
        int nSnapShot;
        int snapInterval;
        int snapTimeBeginItr;

        // DataGrid parameter
        int nData;
        double xOrigin, yOrigin, zOrigin;
        int nxData, nyData, nzData;
        double lxData, lyData, lzData;
        double dxData, dyData, dzData;
        int nNodesInCellData;
        int nCellsDataGlobal;

        // Voxel creation
        int nxOpt, nyOpt, nzOpt;
        double lxOpt, lyOpt, lzOpt;
        double dxOpt, dyOpt, dzOpt;
        int nCellsOptGlobal;
        int nNodesOptGlobal;
        int stepMax;
        std::string inputDir;

        // Boundary parameter for stgrid
        std::vector<std::string> bdStr;
        std::vector<std::string> bdType;
        std::vector<std::vector<double>> bdValue; 

        // For image data
        std::vector<double> phi;

        // For Dirichlet boundary data
        std::vector<int> vDirichletNode;
        std::vector<int> pDirichletNode;
        std::vector<std::vector<double>> vDirichletValue;
        std::vector<double> pDirichletValue;
        std::map<int, std::vector<double>> vDirichlet;
        std::map<int, std::vector<double>> vDirichletWall;
        std::map<int, double> pDirichlet;

        std::vector<int> controlBoundaryMap;
        std::vector<int> planeDir;

        std::vector<bool> isBoundaryEdge;

        // For cell and node data
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> cell;

        std::vector<int> sortCell;
        std::vector<int> sortNode;

        // Error parameter
        bool isReadingError = false;

        // data parameter
        int nControlNodesInCell;
        std::vector<std::vector<std::vector<double>>> velocityData;
        std::vector<std::vector<int>> controlNodeInCell;
        std::vector<int> controlCellMap;

        // post inverse parameter
        int nRef;
        int crossPoint;
        int flowRateVelDir;
        std::vector<std::vector<std::vector<double>>> velRef;
        std::vector<std::vector<std::vector<double>>> velOpt;

        void setSolidBoundary();
        void setFluidDomain();

    private:
        void setApplication(std::string appName);
        void tryOpenConfigFile(std::string inputFile);
        void tryReadConfigFile();
        void readConfigFile();

        void readGridParameter(); 
        void readBoundaryParameter();
        void readControlBoundaryParameter();
        void readStructuredGridParameter();
        void readStructuredBoundaryParameter();
        void readBasicParameter();           
        void readPysicalParameter();
        void readNRParameter();
        void readDarcyParameter();    
        void readTimeParameter();         
        void readInverseParameter();
        void readDataParameter();
        void readPostInverseBasicParameter();
        void readPostInverseVelocityParameter();
        void readPostInverseFlowRateParameter();
        void readVoxelCreationParameter();
        void readBoundaryTypeAndValue(std::string labelType, 
                                      std::string labelValue, int &tmp);
};

#endif