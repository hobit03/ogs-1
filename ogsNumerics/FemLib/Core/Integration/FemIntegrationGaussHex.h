/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationGaussHex.h
 *
 * Created on 2014-03-06 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/Integration/GaussLegendre.h"
#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"

namespace FemLib
{

class FemIntegrationGaussHex : public AbstractFemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        const size_t nGauss = getSamplingLevel();
        size_t gp_r = (igp / (nGauss * nGauss));
        size_t gp_s = (igp % (nGauss * nGauss));
        size_t gp_t = gp_s % nGauss;
        gp_s /= nGauss;
        x[0] = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
        x[1] = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
        x[2] = MathLib::GaussLegendre::getPoint(nGauss, gp_t);
    }

    double getWeight(size_t igp) const
    {
        const size_t nGauss = getSamplingLevel();
        size_t gp_r = (igp / (nGauss * nGauss));
        size_t gp_s = (igp % (nGauss * nGauss));
        size_t gp_t = gp_s % nGauss;
        gp_s /= nGauss;
        return  MathLib::GaussLegendre::getWeight(nGauss, gp_r)
                *MathLib::GaussLegendre::getWeight(nGauss, gp_s)
                *MathLib::GaussLegendre::getWeight(nGauss, gp_t);
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement&, size_t n_sampl_level) const
    {
        return n_sampl_level*n_sampl_level*n_sampl_level;
    }
};


}
