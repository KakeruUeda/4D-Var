#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include "TextParser.h"
#include "Array.h"
#include "MyMPI.h"
extern MyMPI mpi;

enum class Application
{
    STRGRID = 0, SNS = 1, USNS = 2,
    TDVAR = 3,  FDVAR = 4
};

class Config
{
    public:
        Config(std::string inputFile, std::string appName);
        ~Config(){};

        TextParser tp;
        Application app;

        // Basic parameter
        int dim, nOMP;
        std::string outputDir, gridTypeString;

        // Physical parameter
        double rho, mu, Re;

        // Grid parameter
        int nx, ny, nz;
        double lx, ly, lz;
        double dx, dy, dz;
        int nxNodes, nyNodes, nzNodes;
        int nxCells, nyCells, nzCells;
        int nCellsGlobal, nNodesGlobal, nNodesInCell;

        // Boundary parameter for stgrid
        std::vector<std::string> bdStr;
        std::vector<std::string> bdType;
        std::vector<std::vector<double>> bdValue; 

        // For image data
        std::vector<double> phi;

        // For Dirichlet boundary data
        std::vector<int> vDirichletNode;
        std::vector<int> pDirichletNode;
        std::vector<double> vDirichletValue;
        std::vector<double> pDirichletValue;

        // For cell and node data
        std::vector<std::vector<double>> node;
        std::vector<std::vector<double>> cell;

        // Error parameter
        bool isReadingError = false;

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