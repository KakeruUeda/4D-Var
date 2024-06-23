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
#include "define.h"

extern MyMPI mpi;

enum class Application
{
    STRGRID = 0, SNS = 1, USNS = 2,
    TDVAR = 3,  FDVAR = 4
};

enum class GridType
{
    STRUCTURED = 0, 
    UNSTRUCTURED = 1
};

class Config
{
    public:
        Config(std::string inputFile, std::string appName);
        ~Config(){};

        TextParser tp;
        Application app;
        GridType gridType;

        // Basic parameter
        int dim, nOMP;
        std::string outputDir;

        // Physical parameter
        double rho, mu, Re;

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

        // ObservedGrid parameter
        int nxObs, nyObs, nzObs;
        double lxObs, lyObs, lzObs;
        double dxObs, dyObs, dzObs;

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
        std::map<int, double> pDirichlet;

        // For cell and node data
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> cell;

        std::vector<int> sortCell;
        std::vector<int> sortNode;

        // Error parameter
        bool isReadingError = false;

        void setSolidBoundary();
        void setFluidDomain();

    private:
        void setApplication(std::string appName);
        void tryOpenConfigFile(std::string inputFile);
        void tryReadConfigFile();
        void readConfigFile();

        void readGridTypeParameter();  
        void readGridParameter(); 
        void readImageData();
        void readStructuredGridParameter();
        void readStructuredBoundaryParameter();
        void readBasicParameter();           
        void readPysicalParameter(); 
        void readDarcyParameter();    
        void readTimeParameter();         
        void readImageParameter();

        void readBoundaryTypeAndValue(std::string labelType, 
                                      std::string labelValue, int &tmp);
};

#endif