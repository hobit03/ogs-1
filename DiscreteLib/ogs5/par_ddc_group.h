
#pragma once

#include <string>
#include <set>

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "MeshLib/Core/IMesh.h"

#include "par_ddc.h"
#include "rf_num_new.h"
#include "equation_class.h"

namespace OGS5
{

/**
 * \brief a class managing decomposed sub-domains
 */
class CPARDomainGroup
{

public:
    CPARDomainGroup(MeshLib::IMixedOrderMesh &msh, std::set<std::pair<bool, size_t> > &set_property) : _set_property(set_property)
    {
        _msh = &msh;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        for (std::set<std::pair<bool, size_t> >::iterator itr=set_property.begin(); itr!=set_property.end(); ++itr) {
            if (!itr->first) use_linear = true;
            if (itr->first) use_quad = true;
        }
    }

    virtual ~CPARDomainGroup()
    {
        //delete _dom_vector
    }

    void addDomain(CPARDomain* dom) 
    {
        _dom_vector.push_back(dom);
    }

    void setup();

    void DDCAssembleGlobalMatrix();
    void SetBoundaryConditionSubDomain();



private:
    int myrank;
    std::set<std::pair<bool, size_t> > _set_property;
    std::vector<CPARDomain*> _dom_vector;
    std::vector<int> _node_connected_doms;
    MeshLib::IMixedOrderMesh* _msh;
    Linear_EQS* eqs_new;
    CNumerics *num;
    bool use_linear;
    bool use_quad;

    /// Find nodes of all neighbors of each node
    void findNodesOnInterface(bool quadr);
    void assembleGlobalMatrix() {};
    void countDoms2Nodes(bool quad);

    void setupDomain( CPARDomain* m_dom, const std::vector<long>& bc_buffer, bool quadr ) ;

    void solveEQS(size_t problem_id)
    {
        long global_eqs_dim = 1;

        //construct eqs
#if 0
        for (size_t i=0; i<dom_vector.size(); i++) {
            CPARDomain *dom = dom_vector[i];
            dom->InitialEQS(problem_id);
#ifdef USE_MPI
            dom->ConfigEQS(num, global_eqs_dim, msh_order);
#endif
            dom->assembly(problem_id);
        }
#ifndef USE_MPI
        assembleGlobalMatrix();
#endif
#endif
        // solve
#ifdef USE_MPI
        for (size_t i=0; i<_dom_vector.size(); i++) {
            CPARDomain *dom = _dom_vector[i];
            dom->getEQS(false)->Solver(eqs_new->getX(), global_eqs_dim);
        }
#else
        eqs_new->Solver();
#endif

        for (size_t i=0; i<_dom_vector.size(); i++) {
            CPARDomain *dom = _dom_vector[i];
            //dom->CleanEQS();
        }
    }
};

}
