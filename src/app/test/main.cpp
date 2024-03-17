#include "Cell.h"
#include "Config.h"
#include "Direct.h"
#include "PreDirect.h"

int main()
{
    PreDirect pre;

    pre.conf.nx = 2; pre.conf.ny = 4; pre.conf.nz = 6;
    pre.conf.nNodesInCellTmp = 4;
    pre.conf.set();

    Direct direct(pre.conf);

    std::cout << "direct.cell.nCellsGlobal = " << direct.cell.nCellsGlobal << std::endl;
    std::cout << "direct.cell(1).nNodesInCell = " << direct.cell(1).nNodesInCell << std::endl;

}