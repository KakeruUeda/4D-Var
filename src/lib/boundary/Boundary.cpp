#include "Boundary.h"

void StructuredBoundaryFace::setNodesOnBoundaryFace(size_t nxNodes, size_t nyNodes, size_t nzNodes)
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
                                              size_t dim, size_t bdIndex)
{
    dirichletType.resize(node.size());
    dirichletValue.resize(node.size());

    for(int i=0; i<dirichletType.size(); i++)
        dirichletType.at(i) = bdType.at(bdIndex);

    for(size_t i=0; i<dirichletValue.size(); i++)
        for(size_t d=0; d<bdValue.at(bdIndex).size(); d++)
            dirichletValue.at(i).push_back(bdValue.at(bdIndex).at(d));
}


StructuredBoundaryFace createBoundarFaceObject(std::string str)
{
    StructuredBoundaryFace obj(str);
    return obj;
}

//void TopEdgeBoundary::setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz)
//{
//    for(int k=0; k<nz+1; k++)
//        for(int j=0; j<ny+1; j++)
//            for(int i=0; i<nx+1; i++)
//                if(j == ny) this->node.push_back(k * (nx+1) * (ny+1) + j * (nx+1) + i);
//}

//void BottomEdgeBoundary::setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz)
//{
//    for(int k=0; k<nz+1; k++)
//        for(int j=0; j<ny+1; j++)
//            for(int i=0; i<nx+1; i++)
//                if(j == 0) this->node.push_back(k * (nx+1) * (ny+1) + j * (nx+1) + i);
//}

//void LeftEdgeBoundary::setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz)
//{
//   for(int k=0; k<nz+1; k++)
//        for(int j=0; j<ny+1; j++)
//            for(int i=0; i<nx+1; i++)
//                if(i == 0) this->node.push_back(k * (nx+1) * (ny+1) + j * (nx+1) + i);
//}

//void RightEdgeBoundary::setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz)
//{
//   for(int k=0; k<nz+1; k++)
//        for(int j=0; j<ny+1; j++)
//            for(int i=0; i<nx+1; i++)
//                if(i == nx) this->node.push_back(k * (nx+1) * (ny+1) + j * (nx+1) + i);
//}

//void FrontEdgeBoundary::setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz)
//{
//   for(int k=0; k<nz+1; k++)
//        for(int j=0; j<ny+1; j++)
//            for(int i=0; i<nx+1; i++)
//                if(k == 0) this->node.push_back(k * (nx+1) * (ny+1) + j * (nx+1) + i);
//}

//void BackEdgeBoundary::setStructuredEdgeBoundary(size_t nx, size_t ny, size_t nz)
//{
//   for(int k=0; k<nz+1; k++)
//        for(int j=0; j<ny+1; j++)
//            for(int i=0; i<nx+1; i++)
//                if(k == nz) this->node.push_back(k * (nx+1) * (ny+1) + j * (nx+1) + i);
//}


