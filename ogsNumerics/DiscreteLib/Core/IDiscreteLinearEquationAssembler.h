/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteLinearEquationAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IDiscreteVector.h"

namespace MeshLib 
{
class IMesh;
}

namespace MathLib
{
    class ILinearEquation;
}

namespace DiscreteLib
{

/**
 * \brief Interface of discrete system assembler classes
 */
class IDiscreteLinearEquationAssembler
{
public:
    typedef DiscreteLib::IDiscreteVector<double> GlobalVector;

    virtual ~IDiscreteLinearEquationAssembler() {};

    /// assembly
    virtual void assembly( const MeshLib::IMesh &msh, MathLib::ILinearEquation &eqs) = 0;
};

}
