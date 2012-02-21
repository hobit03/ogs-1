
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/DiscreteSystem.h"
#include "NumLib/Solution/IProblem.h"


namespace NumLib
{

/**
 * \brief Interface to all solution algorithm classes
 */
class ISolutionAlgorithm
{
public:
};

/**
 * \brief Abstract class for time-stepping method
 */
class AbstractTimeSteppingAlgorithm : public ISolutionAlgorithm, public ITransientSystem
{
public:
    /// @param tim Time step function
    AbstractTimeSteppingAlgorithm(ITimeStepFunction &tim) : _tim(&tim) {};

    /// get the time step function
    /// @return Time step function
    ITimeStepFunction* getTimeStepFunction() const {return _tim;};

    /// suggest the next time step
    /// @param time_currrent current time step object
    /// @return the next time
    double suggestNext(const TimeStep &time_current)
    {
        return _tim->next(time_current.getTime());
    }

    /// return if the next time step is the same as the given time
    /// @param time the time step object
    /// @bool true if this process is active with the given time
    bool isAwake(const TimeStep &time)
    {
        return (time.getTime()==_tim->next(time.getTime()));
    }

private:
    ITimeStepFunction* _tim;
};


template<   template <class> class T_TIME_ODE_ASSEMBLER, 
            class T_LINEAR_SOLVER, 
            class T_FEM_PROBLEM_AND_ASSEMBLER >
class SingleStepLinearFEM; // : public AbstractTimeSteppingAlgorithm

/**
 * \brief Solution algorithm for linear transient problems using FEM with Euler time stepping method
 *
 * @tparam T_LINEAR_SOLVER     Linear equation solver class
 * @tparam T_USER_FEM_PROBLEM  FEM problem class
 * @tparam T_USER_ASSEMBLY     Local assembly class
 */
template<   template <class> class T_TIME_ODE_ASSEMBLER, 
            class T_LINEAR_SOLVER, 
            template <class> class T_USER_FEM_PROBLEM, 
            class T_USER_ASSEMBLY >
class SingleStepLinearFEM<T_TIME_ODE_ASSEMBLER, T_LINEAR_SOLVER, T_USER_FEM_PROBLEM<T_USER_ASSEMBLY>> : public AbstractTimeSteppingAlgorithm
{
public:
    typedef T_TIME_ODE_ASSEMBLER<T_USER_ASSEMBLY> UserTimeOdeAssembler;
    typedef T_USER_FEM_PROBLEM<T_USER_ASSEMBLY> UserFemProblem;

    ///
    SingleStepLinearFEM(DiscreteSystem &dis, UserFemProblem &problem) 
        : AbstractTimeSteppingAlgorithm(*problem.getTimeSteppingFunction()), 
        _discrete_system(&dis), _problem(&problem), _element_ode_assembler(problem.getElementAssemlby())
    {
        const size_t n_var = problem.getNumberOfVariables();
        // initialize solution vectors
        _u_n.resize(n_var, 0);
        _u_n1.resize(n_var, 0);
        for (size_t i=0; i<n_var; i++) {
            _u_n1[i] = (FemLib::FemNodalFunctionScalar*) problem.getIC(i);
        }
        // create linear equation
        _eqs_id = _discrete_system->addLinearEquation(*new T_LINEAR_SOLVER);
        IDiscreteLinearEquation *linear_eqs = _discrete_system->getLinearEquation(_eqs_id);
        // define dof map
        DofMapManager* dofManager = linear_eqs->getDofMapManger();
        for (size_t i=0; i<n_var; i++) {
            dofManager->addDoF(problem.getField(i)->getNumberOfNodes());
        }
        dofManager->construct();
    };

    /// solve 
    int solveTimeStep(const TimeStep &t_n1)
    {
        const size_t n_var = _problem->getNumberOfVariables();
        for (size_t i=0; i<n_var; i++) {
            if (_u_n[i]!=0) delete _u_n[i];
            _u_n[i] = _u_n1[i];
        }
        return solve(t_n1, _u_n, _u_n1);
    }

    /// get the current solution
    FemLib::FemNodalFunctionScalar* getCurrentSolution(int var_id)
    {
        return _u_n1[var_id];
    }

    /// get the time ode assembler
    UserTimeOdeAssembler* getTimeODEAssembler()
    {
        return &_element_ode_assembler;
    }

private:
    /// 
	int solve(const TimeStep &t_n1, const std::vector<FemLib::FemNodalFunctionScalar*>& u_n, std::vector<FemLib::FemNodalFunctionScalar*>& u_n1) 
    {
        UserFemProblem* pro = _problem;
        const size_t n_var = pro->getNumberOfVariables();
        for (size_t i=0; i<n_var; i++) {
            u_n1[i] = new FemLib::FemNodalFunctionScalar(*u_n[i]);
        }
		
		// collect data
		const double delta_t = t_n1.getTimeStepSize();

        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) 
            pro->getFemDirichletBC(i)->setup();
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) 
            pro->getFemNeumannBC(i)->setup();

		// assembly
        const MeshLib::IMesh *msh = u_n1[0]->getMesh();
        IDiscreteLinearEquation* discreteEqs = _discrete_system->getLinearEquation(_eqs_id);
        std::vector<size_t> list_dof_id(n_var);
        for (size_t i=0; i<n_var; i++)
            list_dof_id[i] = i;
        discreteEqs->construct(ElementBasedAssembler(t_n1, _element_ode_assembler));

        //apply BC1,2
        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) 
            pro->getFemDirichletBC(i)->apply(*discreteEqs->getLinearEquation());
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) 
            pro->getFemNeumannBC(i)->apply(discreteEqs->getRHS());
		
		// solve
		discreteEqs->solve();
        double *x = discreteEqs->getX();

        for (size_t i=0; i<n_var; i++) {
            u_n1[i]->setNodalValues(x); //TODO use DOF
        }

		// check 

        return 0;
	}

private:
    UserFemProblem* _problem;
    UserTimeOdeAssembler _element_ode_assembler;

    DiscreteSystem *_discrete_system;
    size_t _eqs_id;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n1;
};

}
