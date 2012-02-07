
#pragma once

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/Function/Function.h"
#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/Function/FemFunction.h"
#include "IFemBC.h"

namespace FemLib
{
//----------------------------------------------------------

class IDirichletBCMethod
{
public:
//    virtual void apply(int linearEqs, DirichletBC &bc) = 0;
};

class DiagonalizeMethod : public IDirichletBCMethod
{
public:
    void apply( MathLib::ILinearEquations& eqs, const std::vector<size_t> &vec_nodes, const std::vector<double> &vec_values)
    {
        eqs.setKnownX(vec_nodes, vec_values);
    }
};

class ResizeMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class PenaltyMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class LagrangeMultiplier : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};


/**
 * DirichletBC class
 */
template<typename Tval>
class FemDirichletBC : IFemBC
{
public:
    ///
    explicit FemDirichletBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, MathLib::IFunction<Tval, GeoLib::Point> *bc_func, IDirichletBCMethod *method)
    {
        _var = var;
        _geo = geo;
        _bc_func = bc_func;
        _method = method;
    }
    virtual ~FemDirichletBC()
    {
    }

    /// setup B.C.
    void setup()
    {
        const MeshLib::IMesh *msh = _var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // set values
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = msh->getNodeCoordinatesRef(_vec_nodes[i]);
            _vec_values[i] = _bc_func->eval(*x);
        }
    }

    /// apply B.C.
    void apply( MathLib::ILinearEquations& eqs ) 
    {
        DiagonalizeMethod method;
        method.apply(eqs, _vec_nodes, _vec_values);
    }


private:
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::IFunction<Tval, GeoLib::Point> *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<Tval> _vec_values;
    // node id, var id, value
    IDirichletBCMethod *_method;
};


}
