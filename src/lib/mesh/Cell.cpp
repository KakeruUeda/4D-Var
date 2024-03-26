#include "Cell.h"

inline void CellInfo::setArrayZero(size_t n)
{   
    assert(n != 0);
    node.setZero(); nDofs.setZero(); 
    x.setZero(); y.setZero(); z.setZero(); 
    u.setZero(); v.setZero(); w.setZero();
}
