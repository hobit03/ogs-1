/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalDataType.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <Eigen>
#include "MathLib/LinAlg/LinearEquation/EigenDenseLinearEquation.h"

namespace DiscreteLib
{
typedef MathLib::EigenDenseLinearEquation LocalEquation;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;
typedef Eigen::VectorXd LocalVector;

}
