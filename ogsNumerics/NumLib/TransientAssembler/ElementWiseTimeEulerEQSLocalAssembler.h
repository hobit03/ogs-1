/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTimeEulerEQSLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cassert>

#include "MeshLib/Core/IElement.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "IElementWiseTransientLinearEQSLocalAssembler.h"


namespace NumLib
{

/**
 * \brief Euler scheme element assembler for time ODE formulations
 *
 * @tparam  T_USER_ASSEMBLY     User-given assembler
 */
class ElementWiseTimeEulerEQSLocalAssembler : public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    ElementWiseTimeEulerEQSLocalAssembler() : _theta(1.0)
    {
    };

    virtual ~ElementWiseTimeEulerEQSLocalAssembler() {};

    ///
    void setTheta(double v)
    {
        assert(v>=.0 && v<=1.0);
        _theta = v;
    }

    /// assemble a local linear equation for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param eqs            local algebraic equation
    virtual void assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalEquation &eqs);

protected:
    virtual void assembleODE(const TimeStep &time, const MeshLib::IElement &e, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalMatrix &M, LocalMatrix &K, LocalVector &F)  = 0;

private:
    double _theta;
};



} //end
