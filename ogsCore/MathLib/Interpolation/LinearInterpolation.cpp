/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearInterpolation.cpp
 *
 * Created on 2010-09-07 by Thomas Fischer
 */

#include "LinearInterpolation.h"
#include "BaseLib/binarySearch.h"
#include <algorithm>

#include <iostream>
#include <cassert>

namespace MathLib {

LinearInterpolation::LinearInterpolation(const LinearInterpolation &src)
    : _supporting_points (src._supporting_points), _values_at_supp_pnts (src._values_at_supp_pnts)
{}

LinearInterpolation::LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts)
    : _supporting_points(supporting_points), _values_at_supp_pnts(values_at_supp_pnts)
{
    assert (_supporting_points.size() > 1 && _supporting_points.size() == _values_at_supp_pnts.size());
}

LinearInterpolation::LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts, const std::vector<double>& points_to_interpolate, std::vector<double>& values_at_interpol_pnts)
    : _supporting_points (supporting_points), _values_at_supp_pnts (values_at_supp_pnts)
{
    assert (supporting_points.size() > 1 && supporting_points.size() == values_at_supp_pnts.size());

//    std::cout << "LinearInterpolation::LinearInterpolation support_points, values_at_supp_pnts: " << std::endl;
//    for (size_t k(0); k<supporting_points.size(); k++) {
//        std::cout << supporting_points[k] << " " << values_at_supp_pnts[k] << std::endl;
//    }
//    std::cout << std::endl;
    values_at_interpol_pnts.clear();
    for (size_t k(0); k<points_to_interpolate.size(); k++)
        values_at_interpol_pnts.push_back (this->getValue (points_to_interpolate[k]));
}

LinearInterpolation::~LinearInterpolation()
{}

double LinearInterpolation::getValue ( double pnt_to_interpolate )
{
    // search interval that has the point inside
    size_t interval_idx (std::numeric_limits<size_t>::max());
    for (size_t k(1); k<_supporting_points.size() && interval_idx == std::numeric_limits<size_t>::max(); k++) {
        if (_supporting_points[k-1] <= pnt_to_interpolate && pnt_to_interpolate <= _supporting_points[k]) {
            interval_idx = k-1;
        }
    }

    // compute linear interpolation polynom: y = m * x + n
    long double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));
//    double m ((_values_at_supp_pnts[interval_idx] - _values_at_supp_pnts[interval_idx+1]) / (_supporting_points[interval_idx] - _supporting_points[interval_idx+1]));
//    double n (_values_at_supp_pnts[interval_idx+1] - m * _supporting_points[interval_idx+1]);
    long double n (_values_at_supp_pnts[interval_idx] - m * _supporting_points[interval_idx]);

    return m * pnt_to_interpolate + n;
}

} // end MathLib
