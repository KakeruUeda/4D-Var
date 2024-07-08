#include <iostream>
#include <vector>
#include <fstream>

int main()
{
    int nx = 64;
    int ny = 64;
    int nz = 64;

    double lx = 1e0;
    double ly = 1e0; 
    double lz = 1e0;

    int numOfCells = nx * ny * nz;

    std::vector<double> phi(numOfCells, 0e0);

    int tmp = 0;

    /*
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                if(k >= nz/4 && k < nz*3/4 && i >= nx/4 && i < nx*3/4)
                    phi[tmp] = 1.0;
                tmp++;
            }
        }
    }
    */

    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                phi[tmp] = 1.0;
                tmp++;
            }
        }
    }

    tmp = 0;

    std::ofstream image("image" + std::to_string(nx) + "x" 
                                + std::to_string(ny) + "x" 
                                + std::to_string(nz) + ".dat");
    
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                image << phi[tmp] << std::endl;
                tmp++;
            }
        }
    }
    image.close();

}