/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Newton.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"
#include "DiscreteLib/Core/DiscreteSystem.h"

#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * NewtonRaphson
 */
template <class F_R, class F_DX>
class NewtonRaphson : public INonlinearSolver
{
    DiscreteLib::DiscreteSystem* _dis_sys;
    F_R* _f_r;
    F_DX* _f_dx;
    VectorType* _r;
    VectorType* _dx;

    typedef DiscreteLib::DiscreteVector<double> ImplVector;
public:
    NewtonRaphson(DiscreteLib::DiscreteSystem* dis_sys, F_R* f_r, F_DX* f_dx) : _dis_sys(dis_sys), _f_r(f_r), _f_dx(f_dx)
    {
        _r = _dx = 0;
    };

    virtual ~NewtonRaphson()
    {
        _dis_sys->deleteVector(_r);
        _dis_sys->deleteVector(_dx);
    };

    virtual void solve(const VectorType &x_0, VectorType &x_new)
    {
        if (_r==0) {
            _r = _dis_sys->createVector<ImplVector>(x_0.size());
            _dx = _dis_sys->createVector<ImplVector>(x_0.size());
        }
        MathLib::NRCheckConvergence<VectorType, MathLib::NRErrorAbsResMNormOrRelDxMNorm> check(_option.error_tolerance);
        MathLib::NewtonRaphsonMethod nr;
        nr.solve(*_f_r, *_f_dx, x_0, x_new, *_r, *_dx, _option.max_iteration, &check);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(NewtonRaphson);
};



} //end

