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
        Postprocess(){}
        ~Postprocess(){}

        ObservedGrid obs;

        void extractOutletVelocity(DirectProblem &direct);
        void makeObservedData(DirectProblem &direct);
};

#endif