
#pragma once

#include <cstddef>
#include <iostream>
#include <cmath>

#include "MathLib/Parameter/ParameterSet.h"
#include "MathLib/Coupling/ICoupledProblem.h"
#include "MathLib/Coupling/ParameterProblemMappingTable.h"
#include "PartitionedAlgorithm.h"

namespace MathLib
{

/**
 * \brief Abstract class for iterative partitioned methods
 */
template <class T_CONVERGENCE_CHECK>
class AbstractIterativePartitionedMethod : public IPartitionedAlgorithm
{
public:
	///
    AbstractIterativePartitionedMethod() : _max_itr(100), _epsilon(1e-3) {}
    ///
    AbstractIterativePartitionedMethod(double epsilon, size_t max_count) : _max_itr(max_count), _epsilon(epsilon) {}
    ///
    virtual ~AbstractIterativePartitionedMethod() {};

    /// get max. iteration
    size_t getMaximumIterationCounts() const {return _max_itr;};

    /// get epsilon
    double getEpsilon() const {return _epsilon;};

    /// get iteration count
    size_t getIterationCounts() const {return _itr_count;};

    /// solve
    int solve(const std::vector<ICoupledSystem*> &list_coupled_problems, UnnamedParameterSet &parameter_table, const ParameterProblemMappingTable &mapping);

protected:
    virtual void doPostAfterSolve( ICoupledSystem& /*solution*/, UnnamedParameterSet& /*vars*/, const ParameterProblemMappingTable& /*mapping*/ )  {}

    virtual void doPostAfterSolveAll( UnnamedParameterSet& /*vars*/, const ParameterProblemMappingTable& /*mapping*/ ) {}

    virtual bool isFixed() const = 0;

    /// update parameter table from all problems
    void updateParameterTable(const ParameterProblemMappingTable &mapping, bool is_fixed, UnnamedParameterSet &parameter_table)
    {
        const size_t n_parameters = parameter_table.size();
        for (size_t i=0; i<n_parameters; i++) {
        	const ParameterProblemMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            const ICoupledSystem *solution = sysvarid.first;
            if (solution!=0) {
            	//if (is_fixed)
            	//	parameter_table.setFixed(i, true);
                parameter_table.set(i, *solution->getOutput(sysvarid.second));
            }
        }
    }

    /// update parameter table from one problem
    void updateParameterTable(const ICoupledSystem &src_problem, const ParameterProblemMappingTable &mapping, bool is_fixed, UnnamedParameterSet &parameter_table)
    {
        const size_t n_vars = parameter_table.size();
        for (size_t i=0; i<n_vars; i++) {
            const ParameterProblemMappingTable::PairSysVarId &v = mapping._map_paraId2subproblem[i];
            const ICoupledSystem *tmp_problem = v.first;
            if (tmp_problem==&src_problem) {
            	if (is_fixed)
            		parameter_table.setFixed(i, true);
                parameter_table.set(i, *src_problem.getOutput(v.second));
            }
        }
    }

    /// update problem parameter from parameter table
    void setInputParameters( const UnnamedParameterSet &parameter_table, const std::vector<ParameterProblemMappingTable::PairInputVar> &problem_parameters, ICoupledSystem* problem)
    {
        const size_t n_input_parameters = problem_parameters.size();
        for (size_t j=0; j<n_input_parameters; j++) {
            size_t local_var_id = problem_parameters[j].first;
            size_t shared_var_id = problem_parameters[j].second;
            problem->setInput(local_var_id, parameter_table.get(shared_var_id));
        }
    }

    /// is para
    bool isInputParametersUpdated( const UnnamedParameterSet &parameter_table, const std::vector<ParameterProblemMappingTable::PairInputVar> &problem_parameters, ICoupledSystem* problem)
    {
        const size_t n_input_parameters = problem_parameters.size();
        for (size_t j=0; j<n_input_parameters; j++) {
            size_t local_var_id = problem_parameters[j].first;
            size_t shared_var_id = problem_parameters[j].second;
            if (parameter_table.isUpdated(shared_var_id)) return true;
        }
        return false;
    }

private:
    size_t _max_itr;
    size_t _itr_count;
    double _epsilon;
};

template <class T_CONVERGENCE_CHECK>
int AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>::solve (
		const std::vector<ICoupledSystem*> &list_coupled_problems,
		UnnamedParameterSet &parameter_table,
		const ParameterProblemMappingTable &mapping
		)
{
    const size_t n_subproblems = list_coupled_problems.size();

    // initialize variables
    updateParameterTable(mapping, isFixed(), parameter_table);

    // iteration start
    const size_t max_itr = getMaximumIterationCounts();
    bool is_converged = false;
    size_t i_itr = 0;
    double v_diff = .0;
    UnnamedParameterSet prev_parameter_table;
    do {
        //prev_parameter_table.assign(parameter_table);
        parameter_table.move(prev_parameter_table);

        // compute each solution
        size_t count_calculated = 0;
        for (size_t i=0; i<n_subproblems; i++) {
        	ICoupledSystem *problem = list_coupled_problems[i];
            const std::vector<ParameterProblemMappingTable::PairInputVar> &problem_parameters = mapping._list_subproblem_input_source[i];

            // calculate all anyway in the 1st iteration
            //if (i_itr>0 && !isInputParametersUpdated(parameter_table, problem_parameters, problem))
            //    continue;

            setInputParameters(parameter_table, problem_parameters, problem); 
            problem->solve();
            doPostAfterSolve(*problem, parameter_table, mapping);
            count_calculated++;
        }

        if (count_calculated>0) {
            doPostAfterSolveAll(parameter_table, mapping);
            parameter_table.finishUpdate();
        }

        if (n_subproblems>1) {
            // check convergence
        	T_CONVERGENCE_CHECK check;
            is_converged = check.isConverged(prev_parameter_table, parameter_table, getEpsilon(), v_diff);
            std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
            std::cout.precision(3);
            std::cout << v_diff << " ";
        } else {
            is_converged = true;
        }
        ++i_itr;
    } while (!is_converged && i_itr<max_itr);
    std::cout << std::endl;
    std::cout << "iteration count=" << i_itr << ", error=" << v_diff << std::endl;

    _itr_count = i_itr;
    if (i_itr==max_itr) {
        std::cout << "the iteration reached the maximum count " << max_itr << std::endl;
    }

    return 0;
}


} //end
