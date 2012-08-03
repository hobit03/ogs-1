/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file BruteForceClosestPair.cpp
 *
 * Created on 2011-01-25 by Thomas Fischer
 */

#include "BruteForceClosestPair.h"
#include "MathTools.h"

namespace GeoLib {

BruteForceClosestPair::BruteForceClosestPair(
        std::vector<GeoLib::Point*> const & pnts, size_t& id0, size_t& id1) :
    ClosestPair (pnts, id0, id1)
{
    double sqr_shortest_dist (GeoLib::sqrDist (_pnts[0], _pnts[1]));

    const size_t n_pnts (_pnts.size());
    for (size_t i(0); i<n_pnts; i++) {
        for (size_t j(i+1); j<n_pnts; j++) {
            double sqr_dist (GeoLib::sqrDist (_pnts[i], _pnts[j]));
            if (sqr_dist < sqr_shortest_dist) {
                sqr_shortest_dist = sqr_dist;
                _id0 = i;
                _id1 = j;
            }
        }
    }
    id0 = _id0;
    id1 = _id1;
}

} // end namespace GEOLIB
