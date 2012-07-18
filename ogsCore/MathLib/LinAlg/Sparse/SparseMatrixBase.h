#ifndef SPARSEMATRIXBASE_H
#define SPARSEMATRIXBASE_H

#include "../MatrixBase.h"

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE> class SparseMatrixBase : public MatrixBase<IDX_TYPE>
{
public:
    SparseMatrixBase(IDX_TYPE n1, IDX_TYPE n2) : MatrixBase<IDX_TYPE> (n1,n2) {}
    SparseMatrixBase() : MatrixBase<IDX_TYPE> () {}
    virtual void amux(FP_TYPE d, FP_TYPE const * const __restrict__ x, FP_TYPE * __restrict__ y) const = 0;         // y +=d*Ax
    virtual ~SparseMatrixBase() {};
};

} // end namespace MathLib

#endif
