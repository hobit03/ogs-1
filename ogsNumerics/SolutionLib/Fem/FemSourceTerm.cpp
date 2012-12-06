/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemSourceTerm.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemSourceTerm.h"

#include "MeshLib/Tools/Tools.h"


namespace SolutionLib
{

FemSourceTerm::FemSourceTerm(const MeshLib::IMesh *msh, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
: _msh(msh), _geo(geo), _bc_func(func->clone()), _t(.0)
{
    _is_transient = false;
    _do_setup = true;
}

FemSourceTerm::FemSourceTerm(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
: _msh(NULL), _geo(NULL), _bc_func(NULL), _vec_nodes(vec_node_id), _vec_values(vec_node_values), _t(.0)
{
    _is_transient = false;
    _do_setup = false;
}

FemSourceTerm::FemSourceTerm(const FemSourceTerm &src)
: _msh(src._msh), _geo(src._geo), _bc_func(src._bc_func->clone()), _vec_nodes(src._vec_nodes),
  _vec_values(src._vec_values), _t(src._t), _is_transient(src._is_transient),
  _do_setup(src._do_setup)
{
}

FemSourceTerm::~FemSourceTerm()
{
    BaseLib::releaseObject(_bc_func);
}

FemSourceTerm* FemSourceTerm::clone() const
{
    return new FemSourceTerm(*this);
}

void FemSourceTerm::initCurrentTime(double t)
{
    _t = t;
}

/// setup
void FemSourceTerm::setup(size_t order)
{
    if (!_do_setup) return;
    if (!_is_transient) _do_setup = false;

    _msh->setCurrentOrder(order);
    MeshLib::findNodesOnGeometry(_msh, _geo, &_vec_nodes);
    _vec_values.resize(_vec_nodes.size());
    for (size_t i=0; i<_vec_nodes.size(); i++) {
        const GeoLib::Point* x = _msh->getNodeCoordinatesRef(_vec_nodes[i]);
        NumLib::TXPosition pos(_t, x->getData());
        _bc_func->eval(pos, _vec_values[i]);
    }

    if (!_is_transient)
        _do_setup = false;
}

}
