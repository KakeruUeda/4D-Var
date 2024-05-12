#include "DirectProblem.h"

void DirectProblem::prepareSerialMatrix()
{

}

void DirectProblem::prepareParallelMatrix()
{
    grid.divideWholeGrid();
}