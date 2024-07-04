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
        int nControlCellsGlobal, nControlNodesGlobal;
        
        std::vector<double> dirichletBCsValue;
        std::vector<double> dirichletBCsValueNew;
        std::vector<double> dirichletBCsValueInit;
        std::vector<double> dirichletBCsValueNewInit;

        std::vector<std::map<int, std::vector<double>>> vDirichlet;
        std::vector<std::map<int, std::vector<double>>> vDirichletWall;
        std::vector<std::map<int, double>> pDirichlet;
        std::vector<std::map<int, std::vector<double>>> vDirichletNew;
        std::vector<std::map<int, std::vector<double>>> vDirichletWallNew;
        std::vector<std::map<int, double>> pDirichletNew;

        std::vector<int> controlBoundaryMap;
        std::vector<int> controlCellMap;
        std::vector<std::vector<int>> controlNodeInCell;
        
        void initialize(Config &conf);
        void initializeAdjoint(Config &conf);

        void assignDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                                std::vector<std::map<int, double>> &pDirichletNew, Node &node, 
                                int &dim, const int t);
        void assignConstantDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                                        std::vector<std::map<int, double>> &pDirichletNew, Node &node, 
                                        int &dim, const int t);
        void assignPulsatileBCs(const int &t, const double &dt, 
                                const double &T, const int &nDofsGlobal);
        void applyDirichletBCs(Cell &cell, PetscSolver &petsc);
        void applyDirichletBCsAdjoint(Cell &cell, PetscSolver &petsc);
    
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
