#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include <cassert>
#include "Array.h"
#include "Config.h"

template <class T> 
class Node
{
    public:
        Node(){};
        Node(Config &conf) :
        nNodesGlobal(conf.nNodesGlobal), data(conf.nNodesGlobal)
        {  
            for(size_t i=0; i<nNodesGlobal; i++) 
                data[i].u.resize(conf.dim);
            for(size_t i=0; i<nNodesGlobal; i++) 
                data[i].p = 0;
            for(size_t i=0; i<nNodesGlobal; i++) 
                data[i].x.resize(conf.dim);
        }
        ~Node(){}
        
        inline T& operator()(size_t n)
        { return data[n]; }

        inline void resize(size_t n)
        { data.resize(n); }

        inline void setZero()
        { for(T& n : data) n = 0; }

        size_t nNodesGlobal;

    private:
	    std::vector<T> data;
};


struct NodeInfo
{
    public:
        Array1D<size_t> nDofs;
        Array1D<double> x;
        Array1D<double> u;
        double p;
};



#endif