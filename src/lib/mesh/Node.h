#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include <cassert>
#include "Array.h"
#include "Config.h"

struct NodeInfo
{
    public:
        int subId, nDofsOnNode;
        Array1D<int> dofNum, dofNumBd;
        Array1D<double> x, u;
        double p;

        Array1D<int> isDofDirichlet;
};

class Node
{
    public:
        Node(){};
        Node(Config conf) :
        nNodesGlobal(conf.nNodesGlobal), 
        nDofsGlobal(0), data(conf.nNodesGlobal){}
        ~Node(){}
        
        int nNodesGlobal, nDofsGlobal;

        inline NodeInfo& operator()(int n)
        { return data[n]; }

        inline void resize(int n)
        { data.resize(n); }

        void initialize(Config conf);

    private:
	    std::vector<NodeInfo> data;
};


#endif