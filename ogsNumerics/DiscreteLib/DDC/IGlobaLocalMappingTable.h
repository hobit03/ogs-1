/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IGlobaLocalMappingTable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>

namespace DiscreteLib
{

/**
 * \brief Interface of mapping tables between global and local index
 */
class IGlobaLocalMappingTable
{
public:
    virtual ~IGlobaLocalMappingTable() {};
    virtual bool hasGlobal(size_t global_id) = 0;
    virtual size_t global2local(size_t global) = 0;
    virtual size_t local2global(size_t local) = 0;
};

} //end

