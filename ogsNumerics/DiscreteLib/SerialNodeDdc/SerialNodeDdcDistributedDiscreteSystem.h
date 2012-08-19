/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodeDDCSerialSharedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/DiscreteEnums.h"
#include "TemplateSerialNodeDdcDiscreteSystem.h"
#include "SerialNodeDdcDistributedLinearEquation.h"

namespace DiscreteLib
{

/// \brief Discrete system for node-based decomposition using distributed memory
typedef TemplateSerialNodeDdcDiscreteSystem<SerialNodeDdcDistributedLinearEquation, DiscreteSystemType::SerialDistributed> SerialNodeDdcDistributedDiscreteSystem;

} //end
