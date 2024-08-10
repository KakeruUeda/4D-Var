#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

int main()
{
    int nx = 128;
    int ny = 64;
    int nz = 64;

    int lx = 4.0;
    int ly = 2.0;
    int lz = 2.0;

    int nCells = nx * ny * nz;

    std::vector<double> image;

    std::string str;
    std::string imageFile = "../usns/input/unilateralStenosis/image.dat";
    std::ifstream ifsImage(imageFile);

    while (getline(ifsImage, str))
    {
        std::istringstream iss(str);
        getline(iss, str, ' ');
        image.push_back(stod(str));
    }
    ifsImage.close();

    if (nCells != image.size())
    {
        std::cout << nCells << " " << image.size() << std::endl;
        std::cerr << "error" << std::endl;
    }

    std::vector<double> imageDA;

    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                int n = i + j * nx + k * nx * ny;
                if (i >= 32 && i < 96)
                {
                    imageDA.push_back(image[n]);
                }
            }
        }
    }

    std::ofstream out("image.dat");
    for (int ic = 0; ic < imageDA.size(); ic++)
    {
        out << imageDA[ic] << std::endl;
    }
    out.close();
}