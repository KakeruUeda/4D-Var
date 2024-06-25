#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#include <vector>
#include <memory>
#include "Array.h"
#include "Node.h"
#include "Cell.h"
#include "PetscSolver.h"
#include "Config.h"

class DirichletBoundary
{
    public:
        DirichletBoundary(){}
        DirichletBoundary(Config &conf){}
        virtual ~DirichletBoundary(){}

        int nNodesVelocity, nNodesPressure;
        
        std::vector<double> dirichletBCsValue;
        std::vector<double> dirichletBCsValueNew;
        std::vector<double> dirichletBCsValueInit;
        std::vector<double> dirichletBCsValueNewInit;

        std::map<int, std::vector<double>> vDirichlet;
        std::map<int, double> pDirichlet;
        std::map<int, std::vector<double>> vDirichletNew;
        std::map<int, double> pDirichletNew;

        std::vector<int> controlBoundaryMap;
        
        void initialize(Config &conf);
        void initializeAdjoint(Config &conf);

        void assignDirichletBCs(Node &node, int &dim);
        void assignPulsatileBCs(const double &t, const double &dt, 
                                const double &T, const int &nDofsGlobal);
        void applyDirichletBCs(Cell &cell, PetscSolver &petsc);
    
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
