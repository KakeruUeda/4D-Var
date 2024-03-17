#include "Cell.h"

inline void CellProperty::resizeArrayZero(size_t n)
{   
    assert(n != 0);
    node.setZero(n); nDofs.setZero(n); 
    x.setZero(n); y.setZero(n); z.setZero(n); 
    u.setZero(n); v.setZero(n); w.setZero(n);
}
