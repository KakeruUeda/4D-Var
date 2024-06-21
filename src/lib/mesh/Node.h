#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include <cassert>
#include "Array.h"
#include "Config.h"


class Node
{
    public:
        Node(){};
        Node(Config &conf) :
        nNodesGlobal(conf.nNodesGlobal),
        nNodesLocal(0){}
        virtual ~Node(){}
        
        int nNodesGlobal, nNodesLocal;

        std::vector<int> map, mapNew;
        std::vector<int> subId;
        std::vector<int> nDofsOnNode;
        std::vector<std::vector<int>> dofsMap, dofsBCsMap;
        std::vector<std::vector<int>> dofsMapNew, dofsBCsMapNew;
        std::vector<std::vector<bool>> isDirichlet, isDirichletNew;
        std::vector<std::vector<double>> type;
        std::vector<std::vector<double>> x;
        std::vector<std::vector<double>> v, vPrev;
        std::vector<double> p;

        void initialize(Config &conf);
};


#endif