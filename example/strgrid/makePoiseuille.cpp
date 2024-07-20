#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

int main()
{
    int nx = 64;
    int ny = 64;
    int nz = 64;

    double lx = 2e0;
    double ly = 2e0;
    double lz = 2e0;

    double dx = lx / (double)nx;
    double dy = ly / (double)ny;
    double dz = lz / (double)nz;

    std::string controlFace = "left";

    double center[2] = {1e0, 1e0};

    double R = 0.25;
    double pi = 3.14159265358979323846;

    std::vector<int> node;
    std::vector<double> u, v, w;

    for(int k=0; k<nz+1; k++){
        for(int j=0; j<ny+1; j++){
            for(int i=0; i<nx+1; i++){
                int n = i + j * (nx+1) + k * (nx+1) * (ny+1);
                if(controlFace == "left"){
                    if(i == 0){
                        double rz = center[0] - k * dz;
                        double ry = center[1] - j * dy;
                        double r = sqrt(rz*rz + ry*ry);
                        double vel = 2e0 * (1 - ((r * r) / (R * R)))/(pi * R * R);
                        if(r >= R){
                            node.push_back(n);
                            u.push_back(0e0); 
                            v.push_back(0e0); 
                            w.push_back(0e0); 
                        }
                        if(r < R){
                            node.push_back(n);
                            u.push_back(vel);    
                            v.push_back(0e0); 
                            w.push_back(0e0); 
                        }
                    }
                }
            }
        }
    }

    double Q = 0e0;

    for(int ib=0; ib<node.size(); ib++){
        Q += u[ib] * dy * dz;
    }

    std::cout << "Q = " << Q << std::endl;

    std::ofstream out("velocityDirichletPoiseuille.dat");
    for(int ib=0; ib<node.size(); ib++){
        out << node[ib] << " " << u[ib] << " " << v[ib] << " " << w[ib] << std::endl;
    }
    out.close();
}
