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
        std::vector<int> nDofsOnNode, nDofsOnNodeNew;
        std::vector<std::vector<int>> dofsMap, dofsBCsMap;
        std::vector<std::vector<int>> dofsMapNew, dofsBCsMapNew;
        std::vector<std::vector<bool>> isDirichlet, isDirichletNew;
        std::vector<std::vector<double>> type;
        std::vector<std::vector<double>> x;
        std::vector<std::vector<double>> v, vPrev;
        std::vector<std::vector<double>> lambda;
        std::vector<std::vector<std::vector<double>>> vt;
        std::vector<std::vector<std::vector<double>>> lambdat;
        std::vector<std::vector<double>> pt;
        std::vector<double> p;

        void initialize(Config &conf);
        void initializeNew();
        void initializeAdjoint(Config &conf, std::vector<int> &controlBoundaryMap);
};

class SnapShot
{
    public:
        SnapShot(){}
        SnapShot(Config &conf):
        isSnapShot(conf.isSnapShot),
        nSnapShot(conf.nSnapShot),
        snapInterval(conf.snapInterval),
        snapTimeBeginItr(conf.snapTimeBeginItr){}

        int isSnapShot;
        int nSnapShot;
        int snapInterval;
        int snapTimeBeginItr;

        std::vector<std::vector<std::vector<double>>> v;
        void takeSnapShot(std::vector<std::vector<double>> &_v,
                          const int &snapCount, const int &nNodesGlobal, const int &dim);
};


#endif