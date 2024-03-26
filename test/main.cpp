#include "DirectProblem.h"
#include "PreDirect.h"

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    std::string inputFile = argv[1];
    Direct::Preprocess pre(inputFile);
    if(pre.conf.isReadingError) return -1;

    DirectProblem direct(pre.conf);

    //std::cout << "direct.cell.nCellsGlobal = " << direct.cell.nCellsGlobal << std::endl;
    //std::cout << "direct.cell(1).nNodesInCell = " << direct.cell(1).nNodesInCell << std::endl;

    return 0;
}