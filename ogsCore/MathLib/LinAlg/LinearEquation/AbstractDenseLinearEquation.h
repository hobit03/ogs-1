/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractDenseLinearEquation
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <valarray>
#include <vector>

#include "ILinearEquation.h"
#include "MathLib/LinAlg/Dense/Matrix.h"


namespace MathLib
{

/**
 * \brief 
 */
class AbstractDenseLinearEquation : public ILinearEquation
{
public:
    typedef Matrix<double> MatrixType;
    typedef std::valarray<double> VectorType;

    AbstractDenseLinearEquation() : _A(0) {};

    virtual ~AbstractDenseLinearEquation()
    {
        if (_A) delete _A;
    }

    void create(size_t length, RowMajorSparsity* /*sp*/=0)
    {
        resize(length);
    }

    bool isCreated() const { return _A!=0; };


    void resize(size_t length);


    void reset()
    {
        (*_A) = .0;
        _b.resize(_b.size(), .0);
        _x.resize(_x.size(), .0);
    }

    size_t getDimension() const
    {
        return _A->getNRows();
    }

    MatrixType* getA() 
    {
        return _A;
    }

    double getA(size_t rowId, size_t colId)
    {
        return (*_A)(rowId, colId);
    }

    void setA(size_t rowId, size_t colId, double v)
    {
        (*_A)(rowId, colId) = v;
    }

    void addA(size_t rowId, size_t colId, double v)
    {
        (*_A)(rowId, colId) += v;
    }

    double getRHS(size_t rowId)
    {
        return _b[rowId];
    }

    double* getRHS()
    {
        return &_b[0];
    }

    VectorType* getRHSAsStdVec() {return &_b;};

    void setRHS(size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    void addRHS(size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    double* getX()
    {
        return &_x[0];
    }

    VectorType* getXAsStdVec() {return &_x;};


    void setKnownX(size_t row_id, double x);

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        for (size_t i=0; i<vec_id.size(); ++i)
            setKnownX(vec_id[i], vec_x[i]);
    }

    void printout(std::ostream &os=std::cout) const
    {
        os << "not implemented yet." << std::endl;
    }

private:
    MatrixType *_A;
    VectorType _b;
    VectorType _x;

};

}
