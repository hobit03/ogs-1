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
    //pnt_to_interpolate=2101;
	double s1=0.0, s2=0.0, w=0.0;
	double m_second_guess = 0.0;
    // search interval that has the point inside
    size_t interval_idx (std::numeric_limits<size_t>::max());
    for (size_t k(1); k<_supporting_points.size() && interval_idx == std::numeric_limits<size_t>::max(); k++) {
        if (_supporting_points[k-1] <= pnt_to_interpolate && pnt_to_interpolate <= _supporting_points[k]) {
            interval_idx = k-1;
        }
    }
		
	// 1st intervall (see for-loop above)
 //   if ( interval_idx == 0 )
	//{
 //   // compute linear interpolation polynom: y = m_0 * x + n_0
 //   long double m_0 ((_values_at_supp_pnts[1] - _values_at_supp_pnts[0]) / (_supporting_points[1] - _supporting_points[0]));
 //   long double n_0 (_values_at_supp_pnts[0] - m_0 * _supporting_points[0]);
	//s1 = (0.5 * _supporting_points[interval_idx+1] - 0.5 * _supporting_points[interval_idx+3]) / (0.5 * _values_at_supp_pnts[interval_idx+1] - 0.5 * _values_at_supp_pnts[interval_idx+3]);
	//s2 = (0.5 * _supporting_points[interval_idx+0] - 0.5 * _supporting_points[interval_idx+2]) / (0.5 * _values_at_supp_pnts[interval_idx+0] - 0.5 * _values_at_supp_pnts[interval_idx+2]);
	//w= (_values_at_supp_pnts[interval_idx+1] -_values_at_supp_pnts[interval_idx+2])/(_values_at_supp_pnts[interval_idx+2]-_values_at_supp_pnts[interval_idx+3]);
	////if (w > 1.0) w = 1.0; 
 //   //w = 1.0; 
	//m_second_guess = 1/((1. - w) * s1 + w * s2);
	//m_0 = m_second_guess;
	////return m_0 * pnt_to_interpolate + n_0;             
	//}
	// Out of range (see for-loop above)
    if (interval_idx >= _supporting_points.size() )
	{
		//the ascending case
	    if (_supporting_points[0] <= _supporting_points[_supporting_points.size()-1])
		{
			if (pnt_to_interpolate <= _supporting_points[0])
				return _values_at_supp_pnts[0];
			else if (pnt_to_interpolate >= _supporting_points[_supporting_points.size()-1])
				return _values_at_supp_pnts[ _supporting_points.size() -1];
		}
		//the descending case
	    if (_supporting_points[0] > _supporting_points[_supporting_points.size() -1])
		{
			if (pnt_to_interpolate >= _supporting_points[0])
				return _values_at_supp_pnts[0];
			else if (pnt_to_interpolate <= _supporting_points[_supporting_points.size()-1])
				return _values_at_supp_pnts[_supporting_points.size() -1];
		}		
	}
    // compute linear interpolation polynom: y = m * x + n
    double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));
	double n (_values_at_supp_pnts[interval_idx] - m * _supporting_points[interval_idx]);

 //   if ((interval_idx >= 0) && (interval_idx < _supporting_points.size() - 2))
	//{
	//	s1 = (0.5 * _supporting_points[interval_idx+1] - 0.5 * _supporting_points[interval_idx+3]) / (0.5 * _values_at_supp_pnts[interval_idx+1] - 0.5 * _values_at_supp_pnts[interval_idx+3]);
	//	s2 = (0.5 * _supporting_points[interval_idx] - 0.5 * _supporting_points[interval_idx+2]) / (0.5 * _values_at_supp_pnts[interval_idx] - 0.5 * _values_at_supp_pnts[interval_idx+2]);
	//	w= (_values_at_supp_pnts[interval_idx+1] -_values_at_supp_pnts[interval_idx+2])/(_values_at_supp_pnts[interval_idx+2]-_values_at_supp_pnts[interval_idx+3]);
	//	//if (w > 1.0) w = 1.0;
	//	//if (w < 1.0e-03) w = 1.0e-03;
	//	m_second_guess = 1/((1. - w) * s1 + w * s2);
	//	//m = m_second_guess;
	//}

    return m * pnt_to_interpolate + n;
	
	//// search interval that has the point inside
 //   size_t interval_idx (std::numeric_limits<size_t>::max());
 //   for (size_t k(1); k<_supporting_points.size() && interval_idx == std::numeric_limits<size_t>::max(); k++) {
 //       if (_supporting_points[k-1] <= pnt_to_interpolate && pnt_to_interpolate <= _supporting_points[k]) {
 //           interval_idx = k-1;
 //       }
 //   }
 //   // if the index lower than the first entry, 
 //   // get the value on the first entry. 
 //   if ( interval_idx <= 0 )
 //       return _values_at_supp_pnts[0];
 //   // if the index bigger than the last entry, 
 //   // get the value of the last entry. 
 //   else if (interval_idx >= _supporting_points.size() )
 //       return _values_at_supp_pnts[ _supporting_points.size() -1];
	//// compute linear interpolation polynom: y = m * x + n
 //   long double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));
 //   //    double m ((_values_at_supp_pnts[interval_idx] - _values_at_supp_pnts[interval_idx+1]) / (_supporting_points[interval_idx] - _supporting_points[interval_idx+1]));
 //   //    double n (_values_at_supp_pnts[interval_idx+1] - m * _supporting_points[interval_idx+1]);
 //   long double n (_values_at_supp_pnts[interval_idx] - m * _supporting_points[interval_idx]);
	//return m * pnt_to_interpolate + n;


}

// HS added 18.01.2013
double LinearInterpolation::getSlope ( double pnt_to_interpolate, double weight )
{
    double s1=0.0, s2=0.0, w=0.0;
	double m_second_guess = 0.0;
	
	// search interval that has the point inside
    size_t interval_idx (std::numeric_limits<size_t>::max());
    for (size_t k(1); k<_supporting_points.size() && interval_idx == std::numeric_limits<size_t>::max(); k++) {
        if (_supporting_points[k-1] <= pnt_to_interpolate && pnt_to_interpolate <= _supporting_points[k]) {
            interval_idx = k-1;
        }
    }

				
	// 1st intervall (see for-loop above)
 //   if ( interval_idx == 0 )
	//{
	//	// compute linear interpolation polynom: y = m_0 * x + n_0
	//	double m_0 ((_values_at_supp_pnts[1] - _values_at_supp_pnts[0]) / (_supporting_points[1] - _supporting_points[0]));

	//		s1 = (0.5 * _supporting_points[interval_idx+2] - 0.5 * _supporting_points[interval_idx+4]) / (0.5 * _values_at_supp_pnts[interval_idx+2] - 0.5 * _values_at_supp_pnts[interval_idx+4]);
	//		s2 = (0.5 * _supporting_points[interval_idx+1] - 0.5 * _supporting_points[interval_idx+3]) / (0.5 * _values_at_supp_pnts[interval_idx+1] - 0.5 * _values_at_supp_pnts[interval_idx+3]);
	//		w= (weight -_values_at_supp_pnts[interval_idx+3])/(_values_at_supp_pnts[interval_idx+2]-_values_at_supp_pnts[interval_idx+3]);
	//		//if (w > 1.0) w = 1.0; 
	//		m_second_guess = 1/((1. - w) * s1 + w * s2);
	//		m_0 = m_second_guess;

	//	//return m_0;       
	//}
	// Out of range (see for-loop above)
    if (interval_idx >= _supporting_points.size() )
	{
		//the ascending case
	    if (_supporting_points[0] <= _supporting_points[_supporting_points.size()-1])
		{
			if (pnt_to_interpolate <= _supporting_points[0])
			{
				//s1 = (0.5 * _supporting_points[0+2] - 0.5 * _supporting_points[0+4]) / (0.5 * _values_at_supp_pnts[0+2] - 0.5 * _values_at_supp_pnts[0+4]);
				//s2 = (0.5 * _supporting_points[0+1] - 0.5 * _supporting_points[0+3]) / (0.5 * _values_at_supp_pnts[0+1] - 0.5 * _values_at_supp_pnts[0+3]);
				//w= (weight -_values_at_supp_pnts[0+3])/(_values_at_supp_pnts[0+2]-_values_at_supp_pnts[0+3]);
				//if (w > 1.0) w = 1.0; 
				//m_second_guess = 1/((1. - w) * s1 + w * s2);
				m_second_guess = ((_values_at_supp_pnts[1] - _values_at_supp_pnts[0]) / (_supporting_points[1] - _supporting_points[0]));
				return m_second_guess;
			}
			if (pnt_to_interpolate >= _supporting_points[_supporting_points.size()-1])
			{
				//s1 = (0.5 * _supporting_points[(_supporting_points.size()-3)+1] - 0.5 * _supporting_points[(_supporting_points.size()-3)+3]) / (0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)+1] - 0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)+3]);
				//s2 = (0.5 * _supporting_points[(_supporting_points.size()-3)] - 0.5 * _supporting_points[(_supporting_points.size()-3)+2]) / (0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)] - 0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)+2]);
				//w= (weight -_values_at_supp_pnts[(_supporting_points.size()-3)+2])/(_values_at_supp_pnts[(_supporting_points.size()-3)+1]-_values_at_supp_pnts[(_supporting_points.size()-3)+2]);
				//if (w > 1.0) w = 1.0; 
				//m_second_guess = 1/((1. - w) * s1 + w * s2);
				m_second_guess = ((_values_at_supp_pnts[_supporting_points.size()-2] - _values_at_supp_pnts[_supporting_points.size()-1]) / (_supporting_points[_supporting_points.size()-2] - _supporting_points[_supporting_points.size()-1]));
				return m_second_guess;
			}
		}
		//the descending case
	    if (_supporting_points[0] > _supporting_points[_supporting_points.size() -1])
		{
			if (pnt_to_interpolate >= _supporting_points[0])
			{
				//s1 = (0.5 * _supporting_points[0+1] - 0.5 * _supporting_points[0+3]) / (0.5 * _values_at_supp_pnts[0+1] - 0.5 * _values_at_supp_pnts[0+3]);
				//s2 = (0.5 * _supporting_points[0] - 0.5 * _supporting_points[0+2]) / (0.5 * _values_at_supp_pnts[0] - 0.5 * _values_at_supp_pnts[0+2]);
				//w= (weight -_values_at_supp_pnts[0+2])/(_values_at_supp_pnts[0+1]-_values_at_supp_pnts[0+2]);
				//if (w > 1.0) w = 1.0; 
				//m_second_guess = 1/((1. - w) * s1 + w * s2);
				m_second_guess = ((_values_at_supp_pnts[1] - _values_at_supp_pnts[0]) / (_supporting_points[1] - _supporting_points[0]));
				return m_second_guess;
			}
			else if (pnt_to_interpolate <= _supporting_points[_supporting_points.size()-1])
			{
				//s1 = (0.5 * _supporting_points[(_supporting_points.size()-3)+1] - 0.5 * _supporting_points[(_supporting_points.size()-3)+3]) / (0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)+1] - 0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)+3]);
				//s2 = (0.5 * _supporting_points[(_supporting_points.size()-3)] - 0.5 * _supporting_points[(_supporting_points.size()-3)+2]) / (0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)] - 0.5 * _values_at_supp_pnts[(_supporting_points.size()-3)+2]);
				//w= (weight -_values_at_supp_pnts[(_supporting_points.size()-3)+2])/(_values_at_supp_pnts[(_supporting_points.size()-3)+1]-_values_at_supp_pnts[(_supporting_points.size()-3)+2]);
				//if (w > 1.0) w = 1.0; 
				//m_second_guess = 1/((1. - w) * s1 + w * s2);
				m_second_guess = ((_values_at_supp_pnts[1] - _values_at_supp_pnts[0]) / (_supporting_points[1] - _supporting_points[0]));
				return m_second_guess;
			}
		}		
	}

    // compute linear interpolation polynom: y = m * x + n
    double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));	
	
 //   if ((interval_idx >= 0) && (interval_idx < _supporting_points.size() - 2))
	//{
	//	s1 = (0.5 * _supporting_points[interval_idx+1] - 0.5 * _supporting_points[interval_idx+3]) / (0.5 * _values_at_supp_pnts[interval_idx+1] - 0.5 * _values_at_supp_pnts[interval_idx+3]);
	//	s2 = (0.5 * _supporting_points[interval_idx] - 0.5 * _supporting_points[interval_idx+2]) / (0.5 * _values_at_supp_pnts[interval_idx] - 0.5 * _values_at_supp_pnts[interval_idx+2]);
	//	w= (weight -_values_at_supp_pnts[interval_idx+2])/(_values_at_supp_pnts[interval_idx+2]-_values_at_supp_pnts[interval_idx+3]);
	//	//if (w > 1.0) w = 1.0; 
	//	m_second_guess = 1/((1. - w) * s1 + w * s2);
	//	//m = m_second_guess;
	//}
	return m; 


	//// search interval that has the point inside
 //   size_t interval_idx (std::numeric_limits<size_t>::max());
 //   for (size_t k(1); k<_supporting_points.size() && interval_idx == std::numeric_limits<size_t>::max(); k++) {
 //       if (_supporting_points[k-1] <= pnt_to_interpolate && pnt_to_interpolate <= _supporting_points[k]) {
 //           interval_idx = k-1;
 //       }
 //   }

	//// if the index lower than the first entry, 
 //   // return a slop value of zero. 
 //   if ( interval_idx <= 0 )
 //       return 0.0;
 //   // if the index bigger than the last entry, 
 //   // also return a slop value of zero
 //   else if (interval_idx >= _supporting_points.size() )
 //       return 0.0;

 //   // compute linear interpolation polynom: y = m * x + n
 //   double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));

 //   return m; 


}



} // end MathLib
