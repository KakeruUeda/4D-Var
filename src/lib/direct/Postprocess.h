/**
 * @file Postprocess.h
 * @author k.ueda
 * @date Jun, 2024
 */

#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include "DirectProblem.h"
#include "Grid.h"

class Postprocess
{
    public:
        Postprocess(Config &conf):
        voxel(conf){}
        ~Postprocess(){}

        DataGrid voxel;

        void extractOutletVelocity(DirectProblem &direct);
        void createData(DirectProblem &direct);
};

#endif