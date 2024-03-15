#ifdef CELL_H
#define CELL_H

class Cell
{
    public:
        Cell(const ConfigTextParser &params) : 
        nNodesInCell(params.nNodesInCell), nDofsNode(param.nDofsNode);

        ~Cell(){};

        size_t nNodesInCell

        Array1D<size_t> nDofsNode
        Array1D<double> node;
};

#endif