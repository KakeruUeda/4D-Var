#include <bits/stdc++.h>

int main()
{
    int nx, ny, nz;
    double lx, ly, lz;

    nx = 32;
    ny = 32;
    nz = 32;

    std::ofstream image("image" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + ".dat");
    for(int k=0; k<nz; k++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int i=0; i<nx; i++)
            {
                double value = 1e0;
                image << i << " " << j << " " << value << std::endl;
            }
        }
    }
    image.close();

    return 0;
}