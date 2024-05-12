#include "Node.h"

void Node::initialize(Config conf)
{
    // initialize variables
    for(int in=0; in<nNodesGlobal; in++){
        data[in].u.resize(conf.dim);
        data[in].x.resize(conf.dim);
        data[in].nDofsOnNode = data[in].u.size() + 1;
        data[in].isDofDirichlet.resize(data[in].nDofsOnNode);
        nDofsGlobal += data[in].nDofsOnNode;
    }
    for(int in=0; in<nNodesGlobal; in++){
        data[in].p = 0e0;
        for(int d=0; d<conf.dim; d++){
            data[in].u(d) = 0e0;
            data[in].x(d) = 0e0;
        }
        for(int p=0; p<data[in].nDofsOnNode; p++){
            data[in].isDofDirichlet(p) = 0;
        }
    }
}