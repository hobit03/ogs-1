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
 * - Polynomial order (default: linear)
 * - IC
 * - BC1, 2
 */
template <class T_DIS_SYS>
class FemVariable
{
public:
    typedef typename FemLib::FemNodalFunctionScalar<T_DIS_SYS>::type MyNodalFunctionScalar;

    /**
     *
     * @param id        variable id
     * @param name      variable name
     * @param order     polynomial order
     */
    FemVariable(size_t id, const std::string &name, FemLib::PolynomialOrder::type order = FemLib::PolynomialOrder::Linear)
    : _id(id), _name(name), _order(order), _f_ic(NULL)
    {
    }

    ///
    ~FemVariable()
    {
        BaseLib::releaseObject(_f_ic);
        BaseLib::releaseObjectsInStdVector(_map_bc1);
        BaseLib::releaseObjectsInStdVector(_map_bc2);
    }

    //----------------------------------------------------------------------
    size_t getID() const {return _id;};
    const std::string& getName() const { return _name;}
    FemLib::PolynomialOrder::type getOrder() const {return _order;};

    //----------------------------------------------------------------------
    void setIC(MyNodalFunctionScalar* ic) { _f_ic = ic; };
    MyNodalFunctionScalar* getIC() const { return _f_ic; };


    //----------------------------------------------------------------------
    void addDirichletBC(FemDirichletBC* bc)
    {
        bc->setOrder(_order);
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
        bc2->setOrder(_order);
        _map_bc2.push_back(bc2);
    }
    size_t getNumberOfNeumannBC() const {return _map_bc2.size();};
    IFemNeumannBC* getNeumannBC(size_t bc_id) const
    {
        return _map_bc2[bc_id];
    };

private:
    size_t _id;
    std::string _name;
    FemLib::PolynomialOrder::type _order;
    MyNodalFunctionScalar* _f_ic;
    std::vector<FemDirichletBC*> _map_bc1;
    std::vector<IFemNeumannBC*> _map_bc2;
};

} //end

