#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#include <vector>
#include <memory>
#include "Array.h"


class DirichletBoundary
{
    public:
        DirichletBoundary(){};
        virtual ~DirichletBoundary(){};
        
        Array1D<size_t> node;
        Array1D<size_t> value;
        Array1D<std::string> type; 

        void applyBouncaryConditions();
        void assignBoundaryConfitions();
    
    private:
};

class StructuredBoundaryFace
{
    public:
        size_t size;
        StructuredBoundaryFace(std::string face) : bdFaceStr(face) {};
        ~StructuredBoundaryFace(){};

        std::vector<size_t> node;
        std::vector<std::string> dirichletType;
        std::vector<std::vector<double>> dirichletValue;

        size_t getNodeSize()
        { return node.size(); };

        void setSize(size_t n)
        { size = n; };

        std::string bdFaceStr;
        void setNodesOnBoundaryFace(size_t nxNodes, size_t nyNodes, size_t nzNodes);

        void setDirichletInfo(std::vector<std::string> bdType, 
                              std::vector<std::vector<double>> bdValue, 
                              size_t dim, size_t bdIndex);

    private:
};


StructuredBoundaryFace createBoundarFaceObject(std::string str);




//class EdgeBoundary
//{
//    public:
//        EdgeBoundary(){};
//        virtual ~EdgeBoundary(){};
//        
//    public:
//        // Pure virtual 
//        virtual size_t getNodeSize() = 0;
//        virtual size_t getNode(size_t i) = 0;
//        virtual std::string getDirichletType(size_t i) = 0;
//        virtual double getDirichletValue(size_t i, size_t d) = 0;
//        virtual std::vector<std::string> getDirichletTypeItself() = 0;
//        virtual std::vector<std::vector<double>> getDirichletValueItself() = 0;
//        virtual void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz) = 0;
//
//        void setDirichletInfo(std::vector<std::string> bdType, 
//                              std::vector<std::vector<double>> bdValue, 
//                              std::vector<std::string> &dirichletType, 
//                              std::vector<std::vector<double>> &dirichletValue,
//                              size_t dim, size_t bdIndex);
//    
//    private:
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};
//
//class LeftEdgeBoundary : public EdgeBoundary
//{
//    public:
//        // Define edge boundary nodes
//        void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz);
//        
//        // Getter 
//        size_t getNodeSize()
//        { return node.size(); };
//
//        size_t getNode(size_t i)
//        { return node.at(i); };
//
//        std::string getDirichletType(size_t i)
//        { return dirichletType.at(i); };
//        double getDirichletValue(size_t i, size_t d)
//        { return dirichletValue.at(i).at(d); };
//
//        std::vector<std::string> getDirichletTypeItself()
//        { return dirichletType; }
//        std::vector<std::vector<double>> getDirichletValueItself()
//        { return dirichletValue; }
//
//    private: 
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};
//
//class RightEdgeBoundary : public EdgeBoundary
//{
//    public:
//        // Define edge boundary nodes
//        void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz);
//        
//        // Getter 
//        size_t getNodeSize()
//        { return node.size(); };
//
//        size_t getNode(size_t i)
//        { return node.at(i); };
//
//        std::string getDirichletType(size_t i)
//        { return dirichletType.at(i); };
//        double getDirichletValue(size_t i, size_t d)
//        { return dirichletValue.at(i).at(d); };
//
//        std::vector<std::string> getDirichletTypeItself()
//        { return dirichletType; }
//        std::vector<std::vector<double>> getDirichletValueItself()
//        { return dirichletValue; }
//
//    private: 
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};
//
//
//class TopEdgeBoundary : public EdgeBoundary
//{
//    public:
//        // Define edge boundary nodes
//        void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz);
//        
//        // Getter 
//        size_t getNodeSize()
//        { return node.size(); };
//
//        size_t getNode(size_t i)
//        { return node.at(i); };
//
//        std::string getDirichletType(size_t i)
//        { return dirichletType.at(i); };
//        double getDirichletValue(size_t i, size_t d)
//        { return dirichletValue.at(i).at(d); };
//
//        std::vector<std::string> getDirichletTypeItself()
//        { return dirichletType; }
//        std::vector<std::vector<double>> getDirichletValueItself()
//        { return dirichletValue; }
//
//    private: 
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};
//
//class BottomEdgeBoundary : public EdgeBoundary
//{
//    public:
//        // Define edge boundary nodes
//        void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz);
//        
//        // Getter 
//        size_t getNodeSize()
//        { return node.size(); };
//
//        size_t getNode(size_t i)
//        { return node.at(i); };
//
//        std::string getDirichletType(size_t i)
//        { return dirichletType.at(i); };
//        double getDirichletValue(size_t i, size_t d)
//        { return dirichletValue.at(i).at(d); };
//
//        std::vector<std::string> getDirichletTypeItself()
//        { return dirichletType; }
//        std::vector<std::vector<double>> getDirichletValueItself()
//        { return dirichletValue; }
//
//    private: 
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};
//
//class FrontEdgeBoundary : public EdgeBoundary
//{
//    public:
//        // Define edge boundary nodes
//        void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz);
//        
//        // Getter 
//        size_t getNodeSize()
//        { return node.size(); };
//
//        size_t getNode(size_t i)
//        { return node.at(i); };
//
//        std::string getDirichletType(size_t i)
//        { return dirichletType.at(i); };
//        double getDirichletValue(size_t i, size_t d)
//        { return dirichletValue.at(i).at(d); };
//
//        std::vector<std::string> getDirichletTypeItself()
//        { return dirichletType; }
//        std::vector<std::vector<double>> getDirichletValueItself()
//        { return dirichletValue; }
//
//    private: 
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};
//
//class BackEdgeBoundary : public EdgeBoundary
//{
//    public:
//        // Define edge boundary nodes
//        void setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz);
//        
//        // Getter 
//        size_t getNodeSize()
//        { return node.size(); };
//
//        size_t getNode(size_t i)
//        { return node.at(i); };
//
//        std::string getDirichletType(size_t i)
//        { return dirichletType.at(i); };
//
//        
//        double getDirichletValue(size_t i, size_t d)
//        { return dirichletValue.at(i).at(d); };
//
//        std::vector<std::string> getDirichletTypeItself()
//        { return dirichletType; }
//        std::vector<std::vector<double>> getDirichletValueItself()
//        { return dirichletValue; }
//
//    private: 
//        size_t size;
//        std::vector<size_t> node;
//        std::vector<std::string> dirichletType;
//        std::vector<std::vector<double>> dirichletValue;
//};

//std::unique_ptr<EdgeBoundary> createBoundaryObject(std::string str);

#endif
