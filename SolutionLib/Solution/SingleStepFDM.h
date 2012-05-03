
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/LinearEquation/MeshBasedDiscreteLinearEquation.h"
#include "DiscreteLib/Utils/SparsityBuilder.h"
#include "FdmLib/FdmFunction.h"
#include "FdmLib/BoundaryConditions.h"
#include "FdmLib/FDM.h"
#include "NumLib/TransientAssembler/ElementWiseTransientLinearEQSAssembler.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/Tools/TemplateTransientLinearFDMFunction.h"
#include "SolutionLib/Tools/TemplateTransientResidualFDMFunction.h"
#include "SolutionLib/Tools/TemplateTransientDxFDMFunction.h"

#include "AbstractTimeSteppingAlgorithm.h"

namespace SolutionLib
{

/**
 * \brief Solution algorithm for linear transient problems using FEM with single time stepping method
 *
 * @tparam T_USER_FEM_PROBLEM  	FEM problem class
 * @tparam T_LINEAR_SOLVER     	Linear equation solver class
 */
template <
	class T_USER_PROBLEM,
	class T_LINEAR_SOLVER
    >
class SingleStepFDM
	: public AbstractTimeSteppingAlgorithm
{
public:
    typedef T_USER_PROBLEM UserFemProblem;
    typedef TemplateTransientLinearFDMFunction<UserFemProblem, typename UserFemProblem::LinearAssemblerType> UserLinearFemFunction;
    typedef TemplateTransientResidualFDMFunction<UserFemProblem, typename UserFemProblem::ResidualAssemblerType> UserResidualFunction;
    typedef TemplateTransientDxFDMFunction<UserFemProblem, typename UserFemProblem::JacobianAssemblerType> UserDxFunction;
    typedef NumLib::TemplateDiscreteNonlinearSolver<UserLinearFemFunction, UserResidualFunction, UserDxFunction> NonlinearSolverType;
    typedef T_LINEAR_SOLVER LinearSolverType;

    /// constructor
    /// - initialize solution vectors
    /// - set up DoFs and equation index
    /// - prepare linear equation and solver
    /// - prepare linear functions
    SingleStepFDM(DiscreteLib::DiscreteSystem* dis, UserFemProblem* problem)
        : AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()), 
          _problem(problem), _discrete_system(dis)
    {
    	//this->caterogorizeGridPoint();

        const size_t n_var = problem->getNumberOfVariables();
        // create dof map
        for (size_t i=0; i<n_var; i++) {
            _dofManager.addVariableDoFs(dis->getMesh()->getID(), 0, problem->getField(i)->getNumberOfNodes());
        }
        _dofManager.construct();
        // initialize solution vectors
        _u_n1.resize(n_var, 0);
        _u_n1[0] = (FdmLib::FdmFunctionScalar*) problem->getIC(0)->clone();
        const size_t n_dofs = _dofManager.getTotalNumberOfActiveDoFs();
        _vec_n0 = dis->createVector<double>(n_dofs);
        _vec_n1 = dis->createVector<double>(n_dofs);
        _vec_n1_0 = dis->createVector<double>(n_dofs);
        _vec_st = dis->createVector<double>(n_dofs);
        FdmLib::FdmFunctionScalar *f_ic = (FdmLib::FdmFunctionScalar*) problem->getIC(0); //TODO one var
        *_vec_n1 = *f_ic->getNodalValues();
        // create linear equation systems
        _linear_solver = new LinearSolverType();
        _linear_eqs = _discrete_system->createLinearEquation<DiscreteLib::TemplateMeshBasedDiscreteLinearEquation, LinearSolverType, DiscreteLib::SparsityBuilderFromNodeConnectivity>(*_linear_solver, _dofManager);
        // setup functions
        _f_linear = new UserLinearFemFunction(problem, problem->getLinearAssembler(), _linear_eqs);
        _f_r = new UserResidualFunction(problem, &_dofManager, problem->getResidualAssembler());
        _f_dx = new UserDxFunction(problem, problem->getJacobianAssembler(), _linear_eqs);
        // setup nonlinear solver
        _f_nonlinear = new NonlinearSolverType(dis, _f_linear, _f_r, _f_dx);
    };

    ///
    virtual ~SingleStepFDM()
    {
        _discrete_system->deleteLinearEquation(_linear_eqs);
        _discrete_system->deleteVector(_vec_n0);
        _discrete_system->deleteVector(_vec_n1);
        _discrete_system->deleteVector(_vec_n1_0);
        _discrete_system->deleteVector(_vec_st);
        Base::releaseObject(_linear_solver);
        Base::releaseObject(_f_linear);
        Base::releaseObject(_f_r);
        Base::releaseObject(_f_dx);
        Base::releaseObject(_f_nonlinear);
    }

    /// solve 
    int solveTimeStep(const NumLib::TimeStep &t_n1)
    {
    	// time step
        double dt = t_n1.getTime() - AbstractTimeSteppingAlgorithm::getTimeStepFunction()->getPrevious();
        NumLib::TimeStep this_t_n1;
        this_t_n1.assign(t_n1);
        this_t_n1.setTimeStepSize(dt);

        // bc1
        std::vector<size_t> list_bc1_eqs_id;
        std::vector<double> list_bc1_val;
        for (size_t i=0; i<_problem->getNumberOfDirichletBC(); i++) {
            FdmLib::FdmDirichletBC<double> *bc1 = _problem->getFdmDirichletBC(i);
            bc1->setup();
            size_t varId = 0; //TODO var id
            std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc1->getListOfBCValues();
            const size_t msh_id = _discrete_system->getMesh()->getID();
            convertToEqsValues(_dofManager, varId, msh_id, list_bc_nodes, list_bc_values, list_bc1_eqs_id, list_bc1_val);
        }

        // st
        std::vector<size_t> list_st_eqs_id;
        std::vector<double> list_st_val;
        for (size_t i=0; i<_problem->getNumberOfNeumannBC(); i++) {
        	FdmLib::FdmNeumannBC<double, double> *bc2 = _problem->getFdmNeumannBC(i);
        	bc2->setup();
            std::vector<size_t> &list_bc_nodes = bc2->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc2->getListOfBCValues();

            size_t varid = 0; //TODO var id
            const size_t msh_id = _problem->getMesh()->getID();
            convertToEqsValues(_dofManager, varid, msh_id, list_bc_nodes, list_bc_values, list_st_eqs_id, list_st_val);
        }
        (*_vec_st) = .0;
        for (size_t i=0; i<list_st_eqs_id.size(); i++) {
        	(*_vec_st)[list_st_eqs_id[i]] = list_st_val[i];
        }

        // setup functions
        _f_linear->reset(&t_n1, _vec_n0);
        _f_r->reset(&t_n1, _vec_n0, _vec_st);
        _f_dx->reset(&t_n1, _vec_n0);

        // initial guess
        *_vec_n1_0 = *_vec_n0;
        for (size_t i=0; i<list_bc1_eqs_id.size(); i++) {
			(*_vec_n1_0)[list_bc1_eqs_id[i]] = list_bc1_val[i];
		}
        
        // solve
        _f_nonlinear->solve(*_vec_n1_0, *_vec_n1);

        _u_n1[0]->setNodalValues(*_vec_n1);

        return 0;
    }

    /// get the current solution
    FdmLib::FdmFunctionScalar* getCurrentSolution(int var_id)
    {
        return _u_n1[var_id];
    }

    ///
    UserFemProblem* getProblem() {return _problem;};

    /// get a linear solver
    LinearSolverType* getLinearEquationSolver() { return _linear_solver; };
    /// get a nonlinear solver
    NonlinearSolverType* getNonlinearSolver() { return _f_nonlinear;};

    ///
    virtual void accept(const NumLib::TimeStep &t)
    {
        AbstractTimeSteppingAlgorithm::accept(t);
        *_vec_n0 = *_vec_n1; //copy current value to previous value
    };

private:
    DISALLOW_COPY_AND_ASSIGN(SingleStepFDM);

    void convertToEqsValues(const DiscreteLib::DofEquationIdTable &eqs_map, size_t var_id, size_t msh_id, const std::vector<size_t> &list_node_id, const std::vector<double> &list_node_values, std::vector<size_t> &list_eqs_id, std::vector<double> &list_eqs_val)
    {
        for (size_t j=0; j<list_node_id.size(); j++) {
            size_t pt_id = list_node_id[j];
            if (eqs_map.isActiveDoF(var_id, msh_id, pt_id)) {
                size_t eqs_id = eqs_map.mapEqsID(var_id, msh_id, pt_id);
                double bc_val = list_node_values[j];
                list_eqs_id.push_back(eqs_id);
                list_eqs_val.push_back(bc_val);
            }
        }
    }

#if 0
    void caterogorizeGridPoint()
    {
        /// Sort neighbor point index
        size_t idx_buff[5], buff1, buff0;
        FdmLib::NeighborPoint_Type nbp_type;

        MeshLib::IMesh* msh = _discrete_system->getMesh();
        const size_t n_nodes = msh->getNumberOfNodes();
        const size_t n_cells = msh->getNumberOfElements();
        std::vector<bool> cell_is_active(n_cells, false);
        std::vector<bool> pt_is_active(n_nodes, false);
        _list_node_type.resize(n_nodes);
        std::vector<std::vector<size_t> > pt2active_cells;


        size_t num_cell_in_use = 0;
        /// Loop over grid cells
        for(size_t i=0; i<n_cells; i++)
        {
        	MeshLib::IElement* e = msh->getElemenet(i);
        	std::vector<size_t> list_e_node_id;
        	std::vector<GeoLib::Point> list_e_node_pt;
        	e->getNodeIDList(list_e_node_id);
        	msh->getListOfNodeCoordinates(list_e_node_id, list_e_node_pt);

            GeoLib::Surface *boundary;
            if(  boundary->isPntInSfc(list_e_node_pt[0].getData())
               &&boundary->isPntInSfc(list_e_node_pt[1].getData())
               &&boundary->isPntInSfc(list_e_node_pt[2].getData())
               &&boundary->isPntInSfc(list_e_node_pt[3].getData()))
            {
               num_cell_in_use++;
               cell_is_active[i] = true;

               // internal points
               for (size_t j=0; j<list_e_node_id.size(); i++) {
            	   const size_t pt_id = list_e_node_id[j];
                   pt_is_active[pt_id] = true;
                   pt2active_cells[pt_id].push_back(i);
               }
            }
        }

        MeshLib::TopologySequentialNodes2Elements node2conn_ele(msh);

        for (size_t i=0; i<n_nodes; i++)
        {
        	if (!pt_is_active[i]) continue;

        	const size_t n_connected_ele = node2conn_ele.getConnectedElements(i).size();
        	if (n_connected_ele<4) {
                _list_node_type[i] = FdmLib::border;
        	} else {
                _list_node_type[i] = FdmLib::intern;
        	}
        }

    }
#endif

private:
    DiscreteLib::DofEquationIdTable _dofManager;
    UserFemProblem* _problem;
    LinearSolverType* _linear_solver;
    DiscreteLib::DiscreteSystem *_discrete_system;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    std::vector<FdmLib::FdmFunctionScalar*> _u_n1;
    UserLinearFemFunction* _f_linear;
    UserResidualFunction* _f_r;
    UserDxFunction* _f_dx;
    NonlinearSolverType* _f_nonlinear;
    MyFemVector *_vec_n0;
    MyFemVector *_vec_n1_0;
    MyFemVector *_vec_n1;
    MyFemVector *_vec_st;
    std::vector<FdmLib::Point_Type> _list_node_type;
};


}
