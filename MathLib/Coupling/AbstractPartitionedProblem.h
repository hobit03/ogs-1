
#pragma once

#include <vector>
#include <iostream>

#include "MathLib/Parameter/ParameterSet.h"
#include "MathLib/Parameter/NamedIOSystem.h"

#include "ICoupledProblem.h"
#include "ParameterProblemMappingTable.h"
#include "Algorithm/PartitionedAlgorithm.h"

namespace MathLib
{

template <class T_PROBLEM>
class AbstractPartitionedProblem : public NamedDynamicIOSystem<T_PROBLEM>
{
public:
	typedef size_t InternalID;
	typedef NamedDynamicIOSystem<T_PROBLEM> BaseClass;

    ///
	AbstractPartitionedProblem() {};

    ///
    virtual ~AbstractPartitionedProblem() {};

    /// check consistency
    virtual bool check() const
    {
        bool flag = true;

        // check for subproblems
        for (size_t i=0; i<_list_subproblems.size(); i++) {
            if (!_list_subproblems[i]->check())
                flag = false;
        }

        if (!BaseClass::isValid()) return false;

        // check consistency in input and output of subproblems
        for (size_t i=0; i<_list_subproblems.size(); i++) {
            const T_PROBLEM* subproblem = _list_subproblems[i];
            // check input parameters required for the subproblem
            const ParameterProblemMappingTable::ListOfInputVar &vec_registered_input_var = _map._list_subproblem_input_source[i];

            std::vector<size_t> list1(vec_registered_input_var.size());
            for (size_t j=0; j<list1.size(); j++) {
                list1[j] = vec_registered_input_var[j].first;
            }
            std::vector<size_t> list2(subproblem->getNumberOfInputParameters());
            for (size_t j=0; j<list2.size(); j++) {
                list2[j] = j; //subproblem->getInputParameterKey(j);
            }
            sort(list1.begin(), list1.end());
            sort(list2.begin(), list2.end());
            std::vector<size_t> list_diff(list1.size()+list2.size());
            std::vector<size_t>::iterator it = set_symmetric_difference(list1.begin(), list1.end(), list2.begin(), list2.end(), list_diff.begin());

            int n_diff = it-list_diff.begin();
            if (n_diff > 0) {
                for (int k=0; k<n_diff; k++) {
                    std::cout << "*** Error: Inconsistent parameter found in subproblem " << i << " and param " << list_diff[k] << std::endl;
                }

                flag = false;
            }
        }

        return flag;
    }


#if 0
    /// add parameter without giving reference
    virtual void addInputParameter(const std::string& key)
    {
    	InternalID n = NamedDynamicIOSystem<T_PROBLEM>::registerInputParameter(key);
        _map._map_paraId2subproblem.push_back(std::make_pair<T_PROBLEM*,size_t>(0, 0));
        return n;
    }

    /// add parameter and reference
    /// @param name variable name
    /// @param sys problem
    /// @param internal_id parameter id in the sys
    /// @return parameter id
    InternalID addOutputParameter(size_t key, T_PROBLEM &sub_problem, ExternalKey para_key_in_sub_problem)
    {

    	InternalID var_id = NamedDynamicIOSystem<T_PROBLEM>::registerOutputParameter(key);
        // set reference
    	const InternalID para_id = sub_problem.getOutputParameterID(para_key_in_sub_problem);
        const Parameter* v = sub_problem.getOutput(para_id);
        NamedDynamicIOSystem<T_PROBLEM>::setOutput(var_id, v);
        // make a link between this and sub-problem variable
        ParameterProblemMappingTable::PairSysVarId parObj = std::make_pair(&sub_problem, para_id);
        _map._map_paraId2subproblem.push_back(parObj);
        // update a list of sub-problems
        if (find(sub_problem)<0) {
            addSubProblem(sub_problem);
        }

        return var_id;
    }
#endif

    void addProblem(T_PROBLEM &subproblem)
    {
        _list_subproblems.push_back(&subproblem);
        _map._list_subproblem_input_source.resize(_list_subproblems.size());

        //return _list_subproblems.size();
    }

    size_t getNumberOfSubProblems() const {return _list_subproblems.size();};

    T_PROBLEM* getProblem(size_t i) {return _list_subproblems[i];};

    void connectParameters()
    {
    	if (!BaseClass::isValid()) return;

    	const size_t n_in_para = this->getNumberOfInputParameters();
    	const size_t n_out_para = this->getNumberOfOutputParameters();
    	const size_t n_sub_prob = this->getNumberOfSubProblems();

        _map._map_paraId2subproblem.resize(n_in_para+n_out_para);

    	// find origin of output
    	for (size_t i=0; i<n_out_para; i++) {
    		const std::string &p_name = BaseClass::getOutputParameterName(i);
    		const int this_var_id = BaseClass::getOutputParameterID(p_name);
            const size_t internal_id = BaseClass::getInternalIDFromOutputID(this_var_id);
    		for (size_t j=0; j<n_sub_prob; j++) {
    			T_PROBLEM* sub_pr = _list_subproblems[j];
    			if (sub_pr->hasOutputParameter(p_name)) {
    				size_t para_id = sub_pr->getOutputParameterID(p_name);
    		        const Parameter* v = sub_pr->getOutput(para_id);
    		        BaseClass::setOutput(this_var_id, v);
    		        // make a link between this and sub-problem variable
    		        ParameterProblemMappingTable::PairSysVarId parObj = std::make_pair(sub_pr, para_id);
    		        _map._map_paraId2subproblem[internal_id] = parObj;
    				break;
    			}
    		}
    	}

    	// connect input of each sub problems
		for (size_t i=0; i<n_sub_prob; i++) {
			T_PROBLEM* sub_pr = _list_subproblems[i];
			size_t n_sub_in_para = sub_pr->getNumberOfInputParameters();
    		for (size_t j=0; j<n_sub_in_para; j++) {
        		const std::string &p_name = sub_pr->getInputParameterName(j);
        		size_t this_para_id = 0;
        		if (BaseClass::hasInputParameter(p_name)) {
        			this_para_id = BaseClass::getInputParameterID(p_name);
                    this_para_id = BaseClass::getInternalIDFromInputID(this_para_id);
        		} else if (BaseClass::hasOutputParameter(p_name)) {
        			int out_p_id = BaseClass::getOutputParameterID(p_name);
                    this_para_id = BaseClass::getInternalIDFromOutputID(out_p_id);
        		} else {
                    // parameter not found
                    continue;
                }
                _map._list_subproblem_input_source[i].push_back(std::make_pair(j, this_para_id));
    		}
		}
    }

#if 0
    void connectOutput(const std::string& key, T_PROBLEM &subproblem, const std::string& para_key_in_sub_problem)
    {
    }

    /// connect system input and shared variable
    void connectInput(const std::string& key, T_PROBLEM &subproblem, const std::string& para_key_in_sub_problem)
    {
        int subproblem_id = find(subproblem);
        if (subproblem_id>=0) {
        	InternalID this_para_id = 0;
        	if (BaseClass::hasInputParameter(key)) {
        		this_para_id = BaseClass::getInputParameterID(key);
        	} else {
        		this_para_id = BaseClass::getOutputParameterID(key);
        	}
        	InternalID sub_para_id = subproblem.getInputParameterID(para_key_in_sub_problem);
            _map._list_subproblem_input_source[subproblem_id].push_back(std::make_pair(sub_para_id, this_para_id));
        }
    }
#endif

    /// find subproblem
    int find(const T_PROBLEM& sub) const
    {
        typename std::vector<T_PROBLEM*>::const_iterator itr = std::find(_list_subproblems.begin(), _list_subproblems.end(), &sub);
        if (itr!=_list_subproblems.end()) {
            return (itr - _list_subproblems.begin());
        } else {
            return -1;
        }
    }

protected:
    std::vector<T_PROBLEM*> _list_subproblems;
    ParameterProblemMappingTable _map;

};


}
