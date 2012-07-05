
#pragma once

#include <vector>

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"
#include "FemDirichletBC.h"
#include "FemNeumannBC.h"

#include "FemVariable.h"

#include "SolutionLib/DataType.h"

namespace SolutionLib
{

/**
 * \brief Template class for transient linear FEM functions
 *
 * \tparam T_ASSEMBLER
 */
template <
    class T_LOCAL_ASSEMBLER
    >
class TemplateTransientLinearFEMFunction
	: public NumLib::TemplateFunction<SolutionVector, SolutionVector>
{
public:
    typedef T_LOCAL_ASSEMBLER UserLocalAssembler;

    /// constructor
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientLinearFEMFunction(const std::vector<FemVariable*> &list_var, UserLocalAssembler* asssembler, DiscreteLib::IDiscreteLinearEquation* linear_eqs)
        : _local_assembler(asssembler),  _linear_eqs(linear_eqs),
          _t_n1(0), _u_n0(0), _list_var(list_var)
    {
    };

    ///
	virtual ~TemplateTransientLinearFEMFunction() {};

	///
	NumLib::TemplateFunction<SolutionVector,SolutionVector>* clone() const
	{
    	return new TemplateTransientLinearFEMFunction<
    			UserLocalAssembler
    				>(_list_var, _local_assembler, _linear_eqs);
	}

    /// reset property
    void reset(const NumLib::TimeStep* t, SolutionVector* u_n0)
    {
    	this->_t_n1 = const_cast<NumLib::TimeStep*>(t);
    	this->_u_n0 = u_n0;
    };

    /// solve linear equations discretized with FEM
    /// @param u0	initial guess
    /// @param u_n1 new results
    void eval(const SolutionVector &/*u0*/, SolutionVector &u_n1)
    {
    	// input, output
        const NumLib::TimeStep &t_n1 = *this->_t_n1;
        SolutionVector* u_n = this->_u_n0;

        // setup BC
        for (size_t i=0; i<_list_var.size(); i++) {
        	FemVariable* var = _list_var[i];
            for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
                FemDirichletBC* bc1 = var->getDirichletBC(j);
                bc1->setup();
                _linear_eqs->setPrescribedDoF(i, bc1->getListOfBCNodes(), bc1->getListOfBCValues());
            }
            for (size_t j=0; j<var->getNumberOfNeumannBC(); j++) {
                var->getNeumannBC(j)->setup();
            }
        }

        //TODO temporally
        //std::vector<SolutionVector*> vec_un;
        //vec_un.push_back(const_cast<SolutionVector*>(u_n));
        //std::vector<SolutionVector*> vec_un1;
        //vec_un1.push_back(const_cast<SolutionVector*>(&u_n1));

		// assembly
        _linear_eqs->initialize();
        NumLib::ElementWiseTransientLinearEQSAssembler assembler(&t_n1, u_n, &u_n1, _local_assembler);
        _linear_eqs->construct(assembler);

        //apply BC1,2
        for (size_t i=0; i<_list_var.size(); i++) {
        	FemVariable* var = _list_var[i];
            for (size_t j=0; j<var->getNumberOfNeumannBC(); j++) {
                IFemNeumannBC* bc2 = var->getNeumannBC(j);
                _linear_eqs->addRHS(i, bc2->getListOfBCNodes(), bc2->getListOfBCValues(), -1.0);
            }
        }

		// solve
		_linear_eqs->solve();
        _linear_eqs->getX(u_n1);
    }


private:
    UserLocalAssembler *_local_assembler;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    NumLib::TimeStep* _t_n1;
    SolutionVector* _u_n0;
    std::vector<FemVariable*> _list_var;
};


} //end