/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplatePoint.h
 *
 * Created on 2010-01-28 by Thomas Fischer
 */

#ifndef TEMPLATEPOINT_H_
#define TEMPLATEPOINT_H_

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>

#include "GeoObject.h"

namespace GeoLib {

/**
 * \ingroup GEOLIB
 *
 * \brief class-template for points can be instantiated by a numeric type.
 * \param T the coordinate type
 */
template <class T> class TemplatePoint : public GeoObject
{
public:
    /** default constructor */
    TemplatePoint ();

    /** constructor - constructs a TemplatePoint object
     \param x1 value for the first coordinate
     \param x2 value for the second coordinate
     \param x3 value for the third coordinate
     */
    TemplatePoint (T x1, T x2, T x3);

    /** constructor - constructs a TemplatePoint object
     \param x values for three coordinates
     */
    TemplatePoint (T const* x);

    /** virtual destructor */
    virtual ~TemplatePoint() {};

    virtual GEOTYPE getGeoType() const {return GeoLib::POINT;};

    inline TemplatePoint<T> operator- (const TemplatePoint<T> &p) const
    {

        TemplatePoint<T> pt_diff(_x[0]-p._x[0], _x[1]-p._x[1], _x[2]-p._x[2]);
        return pt_diff;
    }

    inline TemplatePoint<T> operator+ (const TemplatePoint<T> &p) const
    {

        TemplatePoint<T> pt_add(_x[0]+p._x[0], _x[1]+p._x[1], _x[2]+p._x[2]);
        return pt_add;
    }

    inline void operator-= (const TemplatePoint<T> &p)
    {
        for (size_t i=0; i<3; i++)
            this->_x[i] -=  p._x[i];
    }

    inline void operator+= (const TemplatePoint<T> &p)
    {
        for (size_t i=0; i<3; i++)
            this->_x[i] +=  p._x[i];
    }

    /** check if the given object equals to this. */
    inline bool operator== (const TemplatePoint<T> &p) const;
/*
    {
        for (size_t i=0; i<3; i++)
            if (_x[i]!=p._x[i]) return false;
        return true;
    }
*/

    /** \brief const access operator
     *  The access to the point coordinates is like the access to a field. Code example:
     * \code
     * Point<double> point (1.0, 2.0, 3.0);
     * double sqrNrm2 = point[0] * point[0] + point[1] * point[1] + point[2] + point[2];
     * \endcode
    */
    const T& operator[] (size_t idx) const {
        assert (idx <= 2);
        return _x[idx];
    }
    /** \brief access operator (see book Effektiv C++ programmieren - subsection 1.3.2 ).
     * \sa const T& operator[] (size_t idx) const
     */
    T& operator[] (size_t idx) {
        return const_cast<T&> (static_cast<const TemplatePoint&> (*this)[idx]);
    }

    /** returns an array containing the coordinates of the point */
    const T* getData () const { return _x; }

    /** returns an array containing the coordinates of the point */
    T* getData () { return _x; }

    /** write point coordinates into stream (used from operator<<)
     * \param os a standard output stream
    */
    virtual void write (std::ostream &os) const {
        os << _x[0] << " " << _x[1] << " " << _x[2] << std::flush;
    }

    /**
     * write point coordinates into string
    */
    virtual std::string write () const {
        std::ostringstream strStream;
        strStream << _x[0] << " " << _x[1] << " " << _x[2];
        return strStream.str();
    }

    /** read point coordinates into stream (used from operator>>) */
    virtual void read (std::istream &is) {
        is >> _x[0] >> _x[1] >> _x[2];
    }

protected:
    T _x[3];
};

template <class T> TemplatePoint<T>::TemplatePoint() :
    GeoObject()
{
    _x[0] = static_cast<T>(0);
    _x[1] = static_cast<T>(0);
    _x[2] = static_cast<T>(0);
}

template <class T> TemplatePoint<T>::TemplatePoint(T x1, T x2, T x3) :
    GeoObject()
{
    _x[0] = x1;
    _x[1] = x2;
    _x[2] = x3;
}

template <class T> TemplatePoint<T>::TemplatePoint (T const* x) :
    GeoObject()
{
    for (size_t k(0); k<3; k++) _x[k] = x[k];
}

/** overload the output operator for class Point */
template <class T>
std::ostream& operator<< (std::ostream &os, const TemplatePoint<T> &p)
{
    p.write (os);
    return os;
}

/** overload the input operator for class Point */
template <class T>
std::istream& operator>> (std::istream &is, TemplatePoint<T> &p)
{
    p.read (is);
    return is;
}

//template <class T>
//bool TemplatePoint<T>::operator== (const TemplatePoint<T> &p) const
//{
//    for (size_t i=0; i<3; i++)
//        if (_x[i]!=p._x[i]) return false;
//    return true;
//}

template <class T>
inline bool TemplatePoint<T>::operator== (const TemplatePoint<T> &p) const
{
    for (size_t i=0; i<3; i++)
        if (_x[i]!=p._x[i]) return false;
    return true;
}

template <>
inline bool TemplatePoint<double>::operator== (const TemplatePoint<double> &p) const
{
    for (size_t i=0; i<3; i++)
        if (fabs(_x[i]-p._x[i])>std::numeric_limits<double>::epsilon()) return false;
    return true;
}

} // end namespace GEO

#endif /* TEMPLATEPOINT_H_ */
