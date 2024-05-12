#include "Grid.h"

Grid::Grid(Config &conf) : 
cell(conf), node(conf), dirichlet(conf), dim(conf.dim),
nNodesGlobal(conf.nNodesGlobal), nCellsGlobal(conf.nCellsGlobal)
{   
    if(conf.gridTypeString == "Structured")
        gridType = GridType::STRUCTURED;
    else if(conf.gridTypeString == "Unstructured")
        gridType = GridType::UNSTRUCTURED;

}

void Grid::setStructuredGrid(const int nxCells, const int nyCells, const int nzCells, 
                             const int nxNodes, const int nyNodes, const int nzNodes, 
                             const double dx, const double dy, const double dz,
                             const int nNodesInCell, const int dim, 
                             Cell &cell, Node &node)
{
    for(int k=0; k<nzCells; k++)
        for(int j=0; j<nyCells; j++)
            for(int i=0; i<nxCells; i++)
                for(int p=0; p<nNodesInCell; p++)
                    cell(k * nxCells * nyCells + j * nxCells + i).node(p) 
                        = structuredGridNodeSet(nxNodes, nyNodes, nzNodes, i, j, k, p);

    for(int k=0; k<nzNodes; k++)
        for(int j=0; j<nyNodes; j++)
            for(int i=0; i<nxNodes; i++)
                for(int d=0; d<dim; d++)
                    node(k * nxNodes * nyNodes + j * nxNodes + i).x(d) 
                        = structuredGridCoordinateSet(dx, dy, dz, i, j, k, d);

    return;
}

int Grid::structuredGridNodeSet(const int nxNodes, const int nyNodes, const int nzNodes,
                                   const int i, const int j, const int k, const int p)
{
    std::vector<int> nodeSet(8);
    nodeSet.at(0) = i   + j*nxNodes     + k*nxNodes*nyNodes;
    nodeSet.at(1) = i+1 + j*nxNodes     + k*nxNodes*nyNodes;
    nodeSet.at(2) = i+1 + (j+1)*nxNodes + k*nxNodes*nyNodes;
    nodeSet.at(3) = i   + (j+1)*nxNodes + k*nxNodes*nyNodes;
    nodeSet.at(4) = i   + j*nxNodes     + (k+1)*nxNodes*nyNodes;
    nodeSet.at(5) = i+1 + j*nxNodes     + (k+1)*nxNodes*nyNodes;
    nodeSet.at(6) = i+1 + (j+1)*nxNodes + (k+1)*nxNodes*nyNodes;
    nodeSet.at(7) = i   + (j+1)*nxNodes + (k+1)*nxNodes*nyNodes;

    return nodeSet.at(p);
}

double Grid::structuredGridCoordinateSet(const double dx, const double dy, const double dz,
                                         const int i, const int j, const int k, const int d)
{
    std::vector<double> coordinateSet(3);
    coordinateSet.at(0) = i * dx;
    coordinateSet.at(1) = j * dy;
    coordinateSet.at(2) = k * dz;

    return coordinateSet.at(d);
}

void Grid::divideWholeGrid()
{
    std::vector<int> cellId;
    std::vector<int> nodeId;

    if(mpi.myId == 0){
        int kk = 0;
        int nparts = mpi.nId;

        int *eptr, *eind;
        eptr = new int[nCellsGlobal+1];
        eind = new int[nCellsGlobal * cell(0).nNodesInCell];

        for(int in=0; in<nCellsGlobal+1; in++)
            eptr[in] = 0;

        for(int ic=0; ic<nCellsGlobal; ic++)
            for(int p=0; p<cell(0).nNodesInCell; p++)
                eind[ic+p] = 0;
    
        for(int ic=0; ic<nCellsGlobal; ic++){
            eptr[ic+1] = (ic+1) * cell(0).nNodesInCell;
            for(int p=0; p<cell(0).nNodesInCell; p++)
                eind[kk+p] = cell(ic).node(p);
            kk += cell(0).nNodesInCell;
        }

        // 8-noded hexa element
        int ncommon_nodes = 4; 

        idx_t objval;
        idx_t options[METIS_NOPTIONS];

        METIS_SetDefaultOptions(options);

        // Specifies the partitioning method.
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;   
             
        // Total communication volume minimization
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;    

        // C-style numbering is assumed that starts from 0.
        options[METIS_OPTION_NUMBERING] = 0;  

        // METIS partition routine
        int ret = METIS_PartMeshDual(&nCellsGlobal, &nNodesGlobal, eptr, eind, NULL, 
                                     NULL, &ncommon_nodes, &nparts, NULL, options, 
                                     &objval, &cellId[0], &nodeId[0]);

        if(ret == METIS_OK)
            std::cout << "\n\n METIS partition routine success "  << std::endl;
        else
            std::cout << " METIS partition routine FAILED "  << std::endl;

        if(eptr) delete[] eptr;
        if(eind) delete[] eind;

        // debug
        std::string vtuFile;
        outputVTU.exportMeshPartitionVTU(vtuFile, node, cell);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&cellId[0], nCellsGlobal, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nodeId[0], nNodesGlobal, MPI_INT, 0, MPI_COMM_WORLD);

    for(int ic=0; ic<nCellsGlobal; ic++)
        cell(ic).subId = cellId[ic];

    for(int in=0; in<nNodesGlobal; in++)
        node(in).subId = nodeId[in];

    MPI_Barrier(MPI_COMM_WORLD);
    nCellsLocal = count(cellId.begin(), cellId.end(), mpi.myId);
    std::cout << " nCellsLocal =  " << nCellsLocal << '\t' << mpi.myId << '\t' << mpi.nId << std::endl;
  
    MPI_Barrier(MPI_COMM_WORLD);

}



