/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GeoObject.h
 *
 * Created on 2010-08-27 by Thomas Fischer
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

namespace GeoLib {

struct GeoObjType
{
    enum type {
        POINT,
        POLYLINE,
        POLYGON,
        SURFACE,
        VOLUME,
        INVALID
    };
};
    

/**
 * \ingroup GEOLIB
 *
 * \brief Base class for classes Point, Polyline, Surface.
 */

class GeoObject {
public:
    GeoObject() {};
    virtual ~GeoObject() {};

    virtual GeoObjType::type getGeoType() const = 0;
};

} // end namespace GEOLIB

#endif /* GEOOBJECT_H_ */
