/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DataType.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/LinAlg/LinearEquations/EigenLinearEquation.h"

namespace NumLib
{
typedef MathLib::EigenDenseLinearEquation LocalEquation;
typedef LocalEquation::LocalMatrix LocalMatrix;
typedef LocalEquation::LocalVector LocalVector;
} //end

