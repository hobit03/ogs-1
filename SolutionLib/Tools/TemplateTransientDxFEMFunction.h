
#pragma once

#include <vector>

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TransientAssembler/ElementWiseTransientDxEQSAssembler.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"


namespace SolutionLib
{

typedef DiscreteLib::DiscreteVector<double> MyFemVector;

/**
 * \brief Template class for transient linear FEM functions
 *
 * \tparam T_FEM_PROBLEM
 * \tparam T_TIME_ODE_ASSEMBLER
 * \tparam T_SPACE_ASSEMBLER
 * \tparam T_LINEAR_SOLVER
 */
template <
	class T_USER_FEM_PROBLEM,
    class T_LOCAL_JACOBIAN_ASSEMBLER
    >
class TemplateTransientDxFEMFunction
	//: public MathLib::TemplateFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_USER_FEM_PROBLEM UserFemProblem;
    typedef T_LOCAL_JACOBIAN_ASSEMBLER UserLocalJacobianAssembler;

    /// constructor
    /// @param problem		Fem problem
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientDxFEMFunction(UserFemProblem* problem, UserLocalJacobianAssembler* asssembler, DiscreteLib::IDiscreteLinearEquation* linear_eqs)
        : _problem(problem), _local_assembler(asssembler),  _linear_eqs(linear_eqs),
          _t_n1(0), _u_n0(0)
    {
    };

    ///
	virtual ~TemplateTransientDxFEMFunction() {};

	///
	NumLib::TemplateFunction<MyFemVector,MyFemVector>* clone() const
	{
    	return new TemplateTransientDxFEMFunction<
    				T_USER_FEM_PROBLEM,
    				T_LOCAL_JACOBIAN_ASSEMBLER
    				>(_problem, _local_assembler, _linear_eqs);
	}

    /// reset property
    void reset(const NumLib::TimeStep* t, MyFemVector* u_n0)
    {
    	this->_t_n1 = const_cast<NumLib::TimeStep*>(t);
    	this->_u_n0 = u_n0;
    };

    /// solve Newton-Raphson equations
    /// @param u_n1	current solution
    /// @param r 	residual
    /// @param du	increment
    void eval(const MyFemVector &u_n1, const MyFemVector &r, MyFemVector &du)
    {
    	// input, output
        const NumLib::TimeStep &t_n1 = *this->_t_n1;
        MyFemVector* u_n = this->_u_n0;

        // prepare data
        UserFemProblem* pro = _problem;

        // setup BC1
        FemVariable* var = pro->getVariable(0);
        for (size_t i=0; i<var->getNumberOfDirichletBC(); i++) {
            FemLib::IFemDirichletBC* bc1 = var->getDirichletBC(i);
            bc1->setup();
            size_t varid = 0; //?
            std::vector<double> bc_value_for_dx(bc1->getListOfBCNodes().size(), .0);
            _linear_eqs->setPrescribedDoF(varid, bc1->getListOfBCNodes(), bc_value_for_dx);
        }

        //TODO temporally
        std::vector<MyFemVector*> vec_un;
        vec_un.push_back(const_cast<MyFemVector*>(u_n));
        std::vector<MyFemVector*> vec_un1;
        vec_un1.push_back(const_cast<MyFemVector*>(&u_n1));

		// assembly
        _linear_eqs->initialize();
        NumLib::ElementWiseTransientDxEQSAssembler assembler(&t_n1, &vec_un, &vec_un1, _local_assembler);
        _linear_eqs->construct(assembler);

        // set residual
        _linear_eqs->addRHS(r, -1.0);

		// solve
		_linear_eqs->solve();
        _linear_eqs->getX(du);
    }


private:
    UserFemProblem* _problem;
    UserLocalJacobianAssembler *_local_assembler;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    NumLib::TimeStep* _t_n1;
    MyFemVector* _u_n0;
};


} //end
