
#pragma once

#include <map>
#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/Function/FemFunction.h"
#include "IFemBC.h"


namespace FemLib
{

/**
 * Neumann BC
 */
template<typename Tval>
class FemNeumannBC : IFemBC
{
public:
    /// 
    FemNeumannBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, MathLib::IFunction<Tval, GeoLib::Point> *func) 
    {
        _var = var;
        _geo = geo;
        _bc_func = func;
    }

    /// setup BC. 
    void setup()
    {
        MeshLib::IMesh *msh = (MeshLib::IMesh*)_var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // distribute to RHS
        _vec_values.resize(_vec_nodes.size());
        if (_var->getDimension()==1) {
            // no need to integrate
            // get discrete values at nodes
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                const GeoLib::Point* x = _vec_nodes[i]->getData();
                _vec_values[i] = _bc_func->eval(*x);
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
                    const GeoLib::Point* x = e->getNode(i_nod)->getData();
                    nodal_val[i_nod] = _bc_func->eval(*x);
                } 
                // compute integrals
                IFiniteElement *fe_edge = _var->getFiniteElement(e);
                std::vector<double> result(edge_nnodes);
                MathLib::Matrix<double> M(edge_nnodes, edge_nnodes);
                M = .0;
                fe_edge->integrateWxN(0, M);
                M.axpy(1.0, &nodal_val[0], 0.0, &result[0]);
                // add into RHS values
                for (size_t k=0; k<edge_nnodes; k++)
                    map_nodeId2val[e->getNodeID(k)] += result[k];
            }
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                _vec_values[i] = map_nodeId2val[_vec_nodes[i]->getNodeID()];
            }
        }
    }

    ///
    template<typename T>
    void apply( T* globalRHS ) 
    {
        for (size_t i=0; i<this->getNumberOfConditions(); i++)
            (*globalRHS)[this->getConditionDoF(i)] += this->getConditionValue(i);
    }


    size_t getNumberOfConditions() const
    {
        return _vec_nodes.size();
    }

    size_t getConditionDoF( size_t i ) const
    {
        return _vec_nodes[i]->getNodeID();
    }

    double getConditionValue( size_t i ) const
    {
        return _vec_values[i];
    }

private:
    // node id, var id, value
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::IFunction<Tval, GeoLib::Point> *_bc_func;
    std::vector<MeshLib::INode*> _vec_nodes;
    std::vector<Tval> _vec_values;
};



}
