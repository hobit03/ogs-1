/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemVariable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "FemLib/Core/PolynomialOrder.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "FemIC.h"
#include "FemDirichletBC.h"
#include "FemNeumannBC.h"

namespace NumLib
{
class ITXFunction;
}

namespace SolutionLib
{

/**
 * \brief FEM variable
 *
 * - var id and name
 * - IC
 * - BC1, 2
 */
class FemVariable
{
public:
    /**
     *
     * @param id        variable id
     * @param name      variable name
     * @param order     polynomial order
     */
    FemVariable(size_t id, const std::string &name, FemLib::PolynomialOrder::type initial_order = FemLib::PolynomialOrder::Linear)
    : _id(id), _name(name), /* _f_ic(NULL)*/ _ic(NULL), _current_order(initial_order)
    {
    }

    ///
    ~FemVariable()
    {
        //BaseLib::releaseObject(_f_ic);
        BaseLib::releaseObject(_ic);
        BaseLib::releaseObjectsInStdVector(_map_bc1);
        BaseLib::releaseObjectsInStdVector(_map_bc2);
    }

    //----------------------------------------------------------------------
    size_t getID() const {return _id;};
    const std::string& getName() const { return _name;}

    //----------------------------------------------------------------------
    void setIC(FemIC* ic) { _ic = ic; };
    FemIC* getIC() const { return _ic; };

    //----------------------------------------------------------------------
    void addDirichletBC(FemDirichletBC* bc)
    {
        _map_bc1.push_back(bc);
    }
    size_t getNumberOfDirichletBC() const {return _map_bc1.size();};
    FemDirichletBC* getDirichletBC(size_t bc_id) const
    {
        return _map_bc1[bc_id];
    };


    //----------------------------------------------------------------------
    void addNeumannBC(IFemNeumannBC* bc2)
    {
        _map_bc2.push_back(bc2);
    }
    size_t getNumberOfNeumannBC() const {return _map_bc2.size();};
    IFemNeumannBC* getNeumannBC(size_t bc_id) const
    {
        return _map_bc2[bc_id];
    };

    void setCurrentOrder(FemLib::PolynomialOrder::type order) {_current_order = order;};
    FemLib::PolynomialOrder::type getCurrentOrder() const {return _current_order;};

private:
    size_t _id;
    std::string _name;
    FemIC* _ic;
    std::vector<FemDirichletBC*> _map_bc1;
    std::vector<IFemNeumannBC*> _map_bc2;
    FemLib::PolynomialOrder::type _current_order;
};

} //end

