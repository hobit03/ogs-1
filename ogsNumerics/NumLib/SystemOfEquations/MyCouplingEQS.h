
#pragma once

#include <map>
#include <valarray>

#include "BaseLib/StringTools.h"
#include "NumLib/Coupling/MonolithicProblem.h"

#include "MyArrayFunction.h"
#include "MyFunction.h"

namespace NumLib
{

class MyCouplingEQS : public AbstractMonolithicSystem<ICoupledSystem>
{
public:
	typedef NumLib::LocalVector ArrayType;
	typedef MyArrayFunction<ArrayType > ParameterType;

	MyCouplingEQS(const std::vector<Variable*> &variables, const std::vector<LinearEquation*> &equations)
	: _active_variables(variables), _equations(equations)
	{
		findInactiveVariables(equations, _inactive_variables);
		_f = new MyFunction(variables, equations);
        AbstractMonolithicSystem<ICoupledSystem>::resizeInputParameter(_inactive_variables.size());
        AbstractMonolithicSystem<ICoupledSystem>::resizeOutputParameter(_active_variables.size());
        for (size_t i=0; i<_inactive_variables.size(); i++)
            AbstractMonolithicSystem<ICoupledSystem>::setInputParameterName(i, _inactive_variables[i]->name);
        for (size_t i=0; i<_active_variables.size(); i++)
            AbstractMonolithicSystem<ICoupledSystem>::setOutputParameterName(i, _active_variables[i]->name);
	};

    void setInitial(std::vector<ArrayType*> &vec_arrays)
    {
        for (size_t i=0; i<_active_variables.size(); i++) {
            ArrayType* arr = vec_arrays[_active_variables[i]->id];
            setOutput(i, new ParameterType(*arr));
        }
    }

    int solve()
    {
        std::map<size_t, ArrayType > inactive_x;
        for (size_t i=0; i<_inactive_variables.size(); i++) {
        	size_t para_id = _inactive_variables[i]->id;
        	const ParameterType *p = getInput<ParameterType>(i);
            inactive_x[para_id] = *p->getArray();
        }

        NumLib::LocalEquation eqs1;
        _f->eval(inactive_x, eqs1);
        eqs1.solve();
        double *u1 = eqs1.getX();
        size_t offset = 0;
        for (size_t i=0; i<_active_variables.size(); i++) {
        	Variable *var = _active_variables[i];
        	//size_t para_id = var->id;
        	ArrayType ret(var->n_dof);
        	for (size_t j=0; j<var->n_dof; j++)
        		ret[j] = u1[offset + j];
        	offset += var->n_dof;

        	setOutput(i, new ParameterType(ret));
        }

        return 0;
    }

    size_t getNumberOfInputParameters() const {return _inactive_variables.size();};
    size_t getNumberOfParameters() const {return _inactive_variables.size() + _active_variables.size();};

private:
    bool isActiveVariable(const Variable* var) const
    {
    	return existVariableInList(_active_variables, var);
    }

    void findInactiveVariables(const std::vector<LinearEquation*> &equations, std::vector<Variable*> &inactive_variables)
    {
    	for (size_t i=0; i<equations.size(); i++) {
    		LinearEquation* eq = equations[i];
    		for (size_t j=0; j<eq->variables.size(); j++) {
    			Variable* var = eq->variables[j];
                if (!isActiveVariable(var)) {
                	if (!existVariableInList(inactive_variables, var))
                		inactive_variables.push_back(var);
                }
    		}
    	}
    }

    bool existVariableInList(const std::vector<Variable*> &list_var, const Variable* var) const
    {
        for (size_t i=0; i<list_var.size(); i++)
            if (list_var[i]->id==var->id) return true;
        return false;
    }

private:
    MyFunction* _f;
	const std::vector<Variable*> _active_variables;
	std::vector<Variable*> _inactive_variables;
	const std::vector<LinearEquation*> _equations;
};

} //end
