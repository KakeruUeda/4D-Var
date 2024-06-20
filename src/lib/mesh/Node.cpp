#include "Node.h"

void Node::initialize(Config &conf)
{
    // initialize variables
    u.resize(nNodesGlobal, std::vector<double>(conf.dim));
    x.resize(nNodesGlobal, std::vector<double>(conf.dim));
    nDofsOnNode.resize(nNodesGlobal);
    isDirichlet.resize(nNodesGlobal);
    isDirichletNew.resize(nNodesGlobal);
    subId.resize(nNodesGlobal);
    map.resize(nNodesGlobal);
    mapNew.resize(nNodesGlobal);
    dofsMap.resize(nNodesGlobal);
    dofsBCsMap.resize(nNodesGlobal);
    dofsMapNew.resize(nNodesGlobal);
    dofsBCsMapNew.resize(nNodesGlobal);

    for(int in=0; in<nNodesGlobal; in++){
        nDofsOnNode[in] = conf.dim + 1;
        isDirichlet[in].resize(nDofsOnNode[in], false);
        isDirichletNew[in].resize(nDofsOnNode[in], false);
    }

    int tmp = 0;

    for(int in=0; in<nNodesGlobal; in++){
        dofsMap[in].resize(nDofsOnNode[in]);
        dofsBCsMap[in].resize(nDofsOnNode[in]);
        dofsMapNew[in].resize(nDofsOnNode[in]);
        dofsBCsMapNew[in].resize(nDofsOnNode[in]);
        for(int id=0; id<nDofsOnNode[in]; id++){
            dofsMap[in][id] = tmp;
            dofsBCsMap[in][id] = tmp;
            dofsMapNew[in][id] = tmp;
            dofsBCsMapNew[in][id] = tmp;
            tmp++;
        }
    }

    for(int in=0; in<nNodesGlobal; in++)
        for(int d=0; d<conf.dim; d++)
            x[in][d] = conf.node[in][d];

    return;
}