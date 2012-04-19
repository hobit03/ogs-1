
#include "FemExtrapolationAverage.h"

#include <algorithm>
#include <cmath>

namespace FemLib
{
// average
void FeExtrapolationAverage::extrapolate(IFiniteElement &fe, const std::vector<MathLib::Vector> &gp_values, std::vector<MathLib::Vector> &nodal_values)
{
    extrapolate<MathLib::Vector>(fe, gp_values, nodal_values);
}

template<typename Tvalue>
void FeExtrapolationAverage::extrapolate(IFiniteElement &fe, const std::vector<Tvalue> &gp_values, std::vector<Tvalue> &nodal_values)
{
    Tvalue ele_avg;
    ele_avg = .0;
    for (size_t j=0; j<gp_values.size(); j++)
        ele_avg += gp_values[j];
    ele_avg /= gp_values.size();

    std::fill(nodal_values.begin(), nodal_values.end(), ele_avg);
}

} //end namespace

