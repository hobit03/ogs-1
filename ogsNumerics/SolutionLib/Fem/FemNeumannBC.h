/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNeumannBC.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"
#include "NumLib/Function/ITXFunction.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"

#include "IFemNeumannBC.h"


namespace SolutionLib
{

/**
 * \brief Neumann BC class for a variable
 *
 */
class FemNeumannBC : public IFemNeumannBC
{
public:
    /// 
    FemNeumannBC(const MeshLib::IMesh *msh, FemLib::LagrangianFeObjectContainer* feObjects, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func);

    ///
    FemNeumannBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values);

    /// clone this object
    FemNeumannBC* clone() const;

    /// setup B.C.
    /// \param order Polynomial order
    virtual void setup(size_t order);

    /// get a list of boundary condition nodes
    virtual std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};

    /// get a list of boundary condition values
    virtual std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    const MeshLib::IMesh* _msh;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    const GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;

};

}
