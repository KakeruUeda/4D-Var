/**
 * @file Cell.h
 * @author K.Ueda
 * @date August, 2024
 */
#ifndef FEM_H
#define FEM_H

#include <iostream>
#include "Array.h"

class FEM
{  
public:
  Array3D<double> v0, v;
  Array1D<double> p;


};

#endif
