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
    STGRID = 0, SNS = 1, USNS = 2,
    TDVAR = 3,  FDVAR = 4
};

class Config
{
    public:
        Config(std::string inputFile, std::string appName);
        ~Config(){};

        TextParser tp;
        Application app;

        // Badsic param
        size_t dim;
        size_t nOMP;
        std::string outputDir;
        std::string gridTypeString;

        // Physical param
        double rho, mu;

        // Grid param
        Array1D<size_t> nx;
        Array1D<double> lx;
        Array1D<double> dx;
        size_t nCellsGlobal;
        size_t nNodesGlobal;
        size_t nNodesInCellTmp;
        Array1D<double> nNodesInCell;

        // Boundary param
        std::vector<std::string> bdStr;
        std::vector<std::string> bdType;
        std::vector<std::vector<size_t>> bdValue; 

        // Image param
        std::vector<double> phi;

        // Error param
        bool isReadingError;

    private:
        void setApplication(std::string appName);
        void tryOpenConfigFile(std::string inputFile);
        void tryReadConfigFile();
        void readConfigFile();

        void readGridTypeParameter();     
        void readBasicParameter();           
        void readPysicalParameter(); 
        void readBoundaryMethodParameter(); 
        void readXFEMParameter();    
        void readDarcyParameter();    
        void readTimeParameter();    
        void readGridParameter();      
        void readBoundaryParameter();    
        void readImageParameter();

        void readBoundaryTypeAndValue(std::string labelType, std::string labelValue, size_t &tmp);
};

#endif