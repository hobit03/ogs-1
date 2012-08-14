/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/DataType.h"
#include "NumLib/Function/IFunction.h"
#include "TXPosition.h"

namespace NumLib
{

/**
 * \brief Interface of any functions in space-time domain
 *
 * This class aims to be an abstract of spatially and temporally distributed data such as
 * physical quantity (e.g. head, stress) and material property (e.g. permeability).
 *
 * TXFunction
 * - is evaluated at particular position in space-time domain
 * - returns scalar or vector value
 * - has some attributes (e.g constant)
 *
 */
class ITXFunction : public IClonable
{
public:
    typedef LocalMatrix DataType;
    
    ///
    ITXFunction() : _is_const(false), _is_temporally_const(false), _is_spatially_const(false) {};

    ///
    virtual ~ITXFunction() {};

    /// evaluate this function at the given position and return vector data

    /// \param x  position in space and time
    /// \param v  evaluated vector
    virtual void eval(const TXPosition /*x*/, DataType &/*v*/) const {};

    /// evaluate this function at the given position and return scalar value.
    ///
    /// Default behavior of this function is to return the 1st component of the vector.
    /// \param x  position in space and time
    /// \param v  evaluated scalar
    virtual void eval(const TXPosition x, double &v) const
    {
        DataType tmp(1,1);
        this->eval(x, tmp);
        v = tmp(0);
    }

    ///
    virtual void eval(const double* x, DataType &v) const
    {
        TXPosition pos(x);
        this->eval(pos, v);
    }

    ///
    virtual void eval(const double* x, double &v) const
    {
        TXPosition pos(x);
        this->eval(pos, v);
    }

    ///
    virtual ITXFunction* clone() const = 0;

    ///
    bool isConst() const {return _is_const;};
    void isConst(bool b) {_is_const = b;};
    ///
    bool isTemporallyConst() const {return _is_temporally_const;};
    void isTemporallyConst(bool b) {_is_temporally_const = b;};
    ///
    bool isSpatiallyConst() const {return _is_spatially_const;};
    void isSpatiallyConst(bool b) {_is_spatially_const = b;};

private:
    bool _is_const;
    bool _is_temporally_const;
    bool _is_spatially_const;
};

/**
 * \brief Constant value
 */
class TXFunctionConstant : public ITXFunction
{
public:
    explicit TXFunctionConstant(double val) : _vec(1,1)
    {
        _vec(0,0) = val;
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(true);
    };
    explicit TXFunctionConstant(const DataType &val) : _vec(val) {};

    virtual ~TXFunctionConstant() {};

    virtual void eval(const TXPosition /*x*/, DataType &v) const { v = _vec;}
    void eval(double &v) const { v = _vec(0,0); };
    void eval(DataType &v) const { v = _vec; };

    virtual TXFunctionConstant* clone() const { return new TXFunctionConstant(_vec); }

    //virtual bool isConst() const {return true;};
    //virtual bool isTemporallyConst() const {return true;};
    //virtual bool isSpatiallyConst() const {return true;};

private:
    DataType _vec;
};


/**
 *
 */
template <typename Tval>
class Multiplication
{
public:
    void doit(const Tval &v1, const Tval &v2, Tval &val)
    {
        val = v1 * v2;
    }
};

/**
 *
 */
template <class T1, class T2, template <typename> class T_OPERATOR>
class TXCompositFunction : public ITXFunction
{
public:
    TXCompositFunction(T1 &f1, T2 &f2) : _f1(&f1), _f2(&f2) {};
    virtual ~TXCompositFunction() {};

    virtual void eval(const TXPosition &x, DataType &val) const
    {
        DataType v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<DataType> op;
        op.doit(v1, v2, val);
    }

    virtual void eval(const double* x, double &val) const
    {
        double v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<double> op;
        op.doit(v1, v2, val);
    }

    virtual TXCompositFunction* clone() const
    {
        return new TXCompositFunction<T1, T2, T_OPERATOR>(*_f1, *_f2);
    }

private:
    T1 *_f1;
    T2 *_f2;
};

template <class F_VECTOR>
class TXVectorFunctionAsColumnData : public ITXFunction
{
public:
    explicit TXVectorFunctionAsColumnData(F_VECTOR* f_vec, size_t column_id)
    : _f_vec(f_vec), _column_id(column_id)
    {
        ITXFunction::isConst(f_vec->isConst());
        ITXFunction::isTemporallyConst(f_vec->isTemporallyConst());
        ITXFunction::isSpatiallyConst(f_vec->isSpatiallyConst());
    };

    virtual ~TXVectorFunctionAsColumnData() {};

    void resetVectorFunction(F_VECTOR* f_vec) {_f_vec = f_vec; };

    virtual void eval(const TXPosition x, DataType &v) const
    {
        DataType tmp_v;
        _f_vec->eval(x, tmp_v);
        if (tmp_v.array().size()>0) {
            v.resize(1,1);
            v(0,0) = tmp_v(_column_id);
        }
    }

    virtual TXVectorFunctionAsColumnData<F_VECTOR>* clone() const
    {
        return new TXVectorFunctionAsColumnData<F_VECTOR>(_f_vec, _column_id);
    }

private:
    F_VECTOR* _f_vec;
    size_t _column_id;
};

} //end


