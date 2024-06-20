#include "Boundary.h"



void StructuredBoundaryFace::setNodesOnBoundaryFace(int nxNodes, int nyNodes, int nzNodes)
{
    if(bdFaceStr == "top"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(j == nyNodes-1){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }

    if(bdFaceStr == "bottom"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(j == 0){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }

    if(bdFaceStr == "left"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(i == 0){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }

    if(bdFaceStr == "right"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(i == nxNodes-1){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }

    if(bdFaceStr == "front"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(k == 0){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }

    if(bdFaceStr == "back"){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(k == nzNodes-1){
                        node.push_back(k * nxNodes * nyNodes + j * nxNodes + i);
                    }
                }
            }
        }
    }
}


void StructuredBoundaryFace::setDirichletInfo(std::vector<std::string> bdType, 
                                              std::vector<std::vector<double>> bdValue, 
                                              int dim, int bdIndex)
{
    dirichletType.resize(node.size());
    dirichletValue.resize(node.size());

    for(int i=0; i<dirichletType.size(); i++)
        dirichletType[i] = bdType[bdIndex];

    for(int i=0; i<dirichletValue.size(); i++)
        for(int d=0; d<bdValue[bdIndex].size(); d++)
            dirichletValue[i].push_back(bdValue[bdIndex][d]);
}

void DirichletBoundary::initialize(Config &conf)
{
    for(int ib=0; ib<conf.vDirichletNode.size(); ib++)
        velocity[ib].node = conf.vDirichletNode[ib];

    for(int ib=0; ib<conf.vDirichletValue.size(); ib++){
        velocity[ib].value.resize(conf.vDirichletValue[ib].size());
        for(int d=0; d<conf.vDirichletValue[ib].size(); d++){
            velocity[ib].value[d] = conf.vDirichletValue[ib][d];
        }
    }

    for(int ib=0; ib<conf.pDirichletNode.size(); ib++)
        pressure[ib].node = conf.pDirichletNode[ib];

    for(int ib=0; ib<conf.pDirichletValue.size(); ib++)
        pressure[ib].value = conf.pDirichletValue[ib];
       
}
