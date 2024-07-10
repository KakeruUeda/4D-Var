/**
 * @file   MyMPI.h
 * @author K.Ueda
 * @date   Jun, 2024
*/

#ifndef MYMPI_H
#define MYMPI_H

#include <iostream>
#include <vector>
#include <mpi.h>

class MyMPI
{
    public:
        MyMPI(){};
        ~MyMPI(){};
        
        int nId, myId;
        
        inline void setSizeAndRank()
        { MPI_Comm_size(MPI_COMM_WORLD, &nId);
          MPI_Comm_rank(MPI_COMM_WORLD, &myId); };

        inline void printSizeAndRank()
        { printf("nId = %d myId = %d \n", nId, myId); }
};


#endif