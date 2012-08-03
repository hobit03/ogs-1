/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IElemenetWiseVectorLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <valarray>

namespace MeshLib
{
class IElement;
}

namespace DiscreteLib
{

/**
 * \brief Interface of all element local assembler classes
 */
template <class T>
class IElemenetWiseVectorLocalAssembler
{
public:
    typedef std::valarray<T> LocalVectorType;

    virtual ~IElemenetWiseVectorLocalAssembler() {};

    /// assemble a local linear equation for the given element
    virtual void assembly(const MeshLib::IElement &e, LocalVectorType &eqs) = 0;
};

} //end
