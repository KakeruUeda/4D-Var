/**
 * @file Postprocess.h
 * @author K.Ueda
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

        void extractOutletVelocity(DirectProblem &direct, std::vector<int> sortNode);
        void createData(DirectProblem &direct);
};

#endif