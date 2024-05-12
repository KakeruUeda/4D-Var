#include <iostream>
#include <vector>
#include <fstream>

int main()
{
    int nx = 4;
    int ny = 4;
    int nz = 4;

    double lx = 1e0;
    double ly = 1e0; 
    double lz = 1e0;

    int numOfCells = nx * ny * nz;

    std::vector<std::vector<double>> phi;
    phi.resize(ny, std::vector<double>(nx, 0e0));

    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                if(k >= nx/4 && k < nx*3/4 && i >= nx/4 && i < nx*3/4)
                    phi[j][i] = 1.0;
            }
        }
    }
    std::ofstream image("image" + std::to_string(nx) + "x" 
                                 + std::to_string(ny) + "x" 
                                 + std::to_string(nz));
    
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                image << phi[j][i] << std::endl;
            }
        }
    }
    image.close();

}