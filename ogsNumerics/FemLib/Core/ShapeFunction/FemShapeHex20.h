/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeHex20.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeHex20 : public TemplateShapeFunction<3, 20> 
{
    virtual void computeShapeFunction(const double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(const double* pt, double* dN) = 0;
};

}
