#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#include <vector>
#include <memory>
#include "Array.h"
#include "Config.h"

struct boundaryInfo
{
    Array1D<int> node;
    Array1D<double> value;
};

class DirichletBoundary
{
    public:
        DirichletBoundary(){}
        DirichletBoundary(Config &conf):
        velocity(conf.vDirichletValue.size()), 
        pressure(conf.pDirichletValue.size()){}
        virtual ~DirichletBoundary(){}
        
        Array1D<boundaryInfo> velocity;
        Array1D<boundaryInfo> pressure;

        void initialize(Config conf);

        void applyBoundaryConditions();
        void assignBoundaryConfitions();
    
    private:
};

class StructuredBoundaryFace
{
    public:
        int size;
        StructuredBoundaryFace(std::string face) : bdFaceStr(face) {};
        ~StructuredBoundaryFace(){};

        std::vector<int> node;
        std::vector<std::string> dirichletType;
        std::vector<std::vector<double>> dirichletValue;

        int getNodeSize()
        { return node.size(); };

        void setSize(int n)
        { size = n; };

        std::string bdFaceStr;
        void setNodesOnBoundaryFace(int nxNodes, int nyNodes, int nzNodes);

        void setDirichletInfo(std::vector<std::string> bdType, 
                              std::vector<std::vector<double>> bdValue, 
                              int dim, int bdIndex);

    private:
};

#endif
