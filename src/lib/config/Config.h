#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <vector>
#include <cassert>

#include "TextParser.h"
#include "Array.h"

enum class GridType
{
    STRUCTURED = 0,
    UNSTRUCTURED = 1
};
 
class Config
{
    public:
        Config(std::string inputFile);
        ~Config(){};

        TextParser tp;
        GridType gridType;

        size_t dim;
        size_t nOMP;
        std::string outputDir;

        double rho, mu;

        Array1D<size_t> nx;
        Array1D<double> lx;
        Array1D<double> dx;

        size_t nCellsGlobal;
        size_t nNodesGlobal;
        size_t nNodesInCellTmp;
        Array1D<double> nNodesInCell;

        bool isReadingError;

    private:
        void tryOpenConfigFile(std::string inputFile);
        void tryReadConfigFile();
        void readConfigFile();

        void readGridType();     void readBase();           
        void readPysicalParam(); void readBoundaryMethod(); 
        void readXFEMParam();    void readDarcyParam();    
        void readTimeParam();    void readGrid();      
        void readBoundary();    void readImage();
};

#endif