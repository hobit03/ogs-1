/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemDirichletBC.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/ITXFunction.h"

namespace SolutionLib
{

/**
 * \brief DirichletBC data for FEM
 *
 * - mesh
 * - geometry
 * - BC function
 */
class FemDirichletBC
{
public:
    ///
    FemDirichletBC(const MeshLib::IMesh* msh, const GeoLib::GeoObject* geo, NumLib::ITXFunction* bc_func)
    {
        _msh = msh;
        _geo = geo;
        _bc_func = bc_func;
        _is_transient = !bc_func->isTemporallyConst();
        _do_setup = true;
    }

    ///
    virtual ~FemDirichletBC()
    {
    }

    /// setup B.C.
    void setup();

    ///
    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};

    ///
    std::vector<double>& getListOfBCValues() {return _vec_values;};

    ///
    bool isTransient() const {return _is_transient;};

private:
    const MeshLib::IMesh* _msh;
    const GeoLib::GeoObject* _geo;
    NumLib::ITXFunction* _bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;
};


}
