
#pragma once

#include <vector>

#include "MathLib/Function/Function.h"
#include "GeoLib/Core/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"
#include "FdmFunction.h"


namespace FdmLib
{

class IFdmBC {};

/**
 * DirichletBC class
 */
template<typename Tval>
class FdmDirichletBC : IFdmBC, public MathLib::SpatialFunction<Tval>
{
public:
    ///
    FdmDirichletBC(TemplateFDMFunction<Tval> *var, GeoLib::GeoObject *geo, bool is_transient, MathLib::SpatialFunction<Tval> *bc_func)
    {
        _var = var;
        _geo = geo;
        _bc_func = (MathLib::SpatialFunction<Tval>*)bc_func->clone();
        _is_transient = is_transient;
        _do_setup = true;
    }

    virtual ~FdmDirichletBC()
    {
        delete _bc_func;
    }

    /// setup B.C.
    void setup()
    {
        if (!_do_setup) return;
        if (!_is_transient) _do_setup = false;

        const MeshLib::IMesh *msh = _var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // set values
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = msh->getNodeCoordinatesRef(_vec_nodes[i]);
            const MathLib::SpatialPosition &tmp = (const MathLib::SpatialPosition) x->getData();
            _bc_func->eval(tmp, _vec_values[i]);
        }
        if (!_is_transient)
            _do_setup = false;
    }

//    /// apply B.C.
//    void apply( MathLib::ILinearEquations& eqs )
//    {
//        DiagonalizeMethod method;
//        method.apply(eqs, _vec_nodes, _vec_values);
//    }

    virtual void eval(const MathLib::SpatialPosition &x, Tval &v)
    {
        _bc_func->eval(x, v);
    }

    MathLib::SpatialFunction<Tval>* clone() const
    {
    	FdmDirichletBC<Tval> *f = new FdmDirichletBC<Tval>(_var, _geo, _is_transient, _bc_func);
        return f;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<Tval>& getListOfBCValues() {return _vec_values;};

private:
    TemplateFDMFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::SpatialFunction<Tval> *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<Tval> _vec_values;
    bool _is_transient;
    bool _do_setup;
};


template<typename Tval, typename Tflux>
class FdmNeumannBC : IFdmBC, public MathLib::SpatialFunction<Tflux>
{
public:
    ///
	FdmNeumannBC(TemplateFDMFunction<Tval> *var, GeoLib::GeoObject *geo, bool is_transient, MathLib::SpatialFunction<Tflux> *func)
    {
        _var = var;
        _geo = geo;
        _bc_func = (MathLib::SpatialFunction<Tflux>*)func->clone();
        _is_transient = is_transient;
        _do_setup = true;
    }

    /// setup BC.
    void setup()
    {
        if (!_do_setup) return;
        if (!_is_transient) _do_setup = false;

        MeshLib::IMesh *msh = (MeshLib::IMesh*)_var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // distribute to RHS
        _vec_values.resize(_vec_nodes.size());
        if (_var->getDimension()==1) {
            // no need to integrate
            // get discrete values at nodes
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                const GeoLib::Point* x = msh->getNodeCoordinatesRef(_vec_nodes[i]);
                _bc_func->eval((const MathLib::SpatialPosition)x->getData(), _vec_values[i]);
            }
        } else {
            // find edge elements on the geo
            std::vector<MeshLib::IElement*> vec_edge_eles;
            MeshLib::findBoundaryElementsOnGeometry(msh, _geo, &vec_edge_eles);
            // for each edge elements found
            std::map<size_t, Tval> map_nodeId2val;
            for (size_t i=0; i<vec_edge_eles.size(); i++) {
                MeshLib::IElement *e = vec_edge_eles[i];
                const size_t edge_nnodes = e->getNumberOfNodes();
                // set values at nodes
                std::vector<double> nodal_val(edge_nnodes);
                for (size_t i_nod=0; i_nod<edge_nnodes; i_nod++) {
                    const GeoLib::Point* x = msh->getNodeCoordinatesRef(e->getNodeID(i_nod));
                    _bc_func->eval((const MathLib::SpatialPosition)x->getData(), nodal_val[i_nod]);
                }
                //TODO compute integrals
                // add into RHS values
//                for (size_t k=0; k<edge_nnodes; k++)
//                    map_nodeId2val[e->getNodeID(k)] += result[k];
            }
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                _vec_values[i] = map_nodeId2val[_vec_nodes[i]];
            }
        }
    }

    ///
    template<typename T>
    void apply( T* globalRHS )
    {
        for (size_t i=0; i<this->getNumberOfConditions(); i++)
            globalRHS[this->getConditionDoF(i)] -= this->getConditionValue(i);
    }


    size_t getNumberOfConditions() const
    {
        return _vec_nodes.size();
    }

    size_t getConditionDoF( size_t i ) const
    {
        return _vec_nodes[i];
    }

    double getConditionValue( size_t i ) const
    {
        return _vec_values[i];
    }

    void eval(const MathLib::SpatialPosition& x, Tflux &v)
    {
        _bc_func->eval(x, v);
    }

    MathLib::SpatialFunction<Tflux>* clone() const
    {
    	FdmNeumannBC<Tval, Tflux> *f = new FdmNeumannBC<Tval, Tflux>(_var, _geo, _is_transient, _bc_func);
        return f;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<Tval>& getListOfBCValues() {return _vec_values;};

private:
    // node id, var id, value
    TemplateFDMFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::SpatialFunction<Tflux> *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<Tval> _vec_values;
    bool _is_transient;
    bool _do_setup;

};

} //end
