/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 2012-09-06 by Haibing Shao
 */

#include "logog.hpp"

#include "MathLib/DataType.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Function/TXFunctionDirect.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "NestedOdeNRIterationStepInitializer.h"
#include "MathLib/ODE/RungeKutta4.h"

template <class T1, class T2>
bool FunctionConcentrations<T1,T2>::initialize(const BaseLib::Options &option)
{
	size_t i;  // index
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

	// get the transformation class instance here
	this->_ReductionKin = femData->m_KinReductScheme; 
	// make sure the reduction scheme is already initialized. 
	if ( !(this->_ReductionKin->IsInitialized()) ) 
	{
		// error msg
	    ERR("While initialize the Global Implicit Reactive Transport Process, the reduction class has not been correctly initialized! ");
		// then stop the program
		exit(1);
	}

	// first get the number of components
	size_t n_Comp = femData->map_ChemComp.size(); 

	// set concentrations of all components as output
	for ( i=0; i < n_Comp; i++ )
		this->setOutputParameterName( i, femData->map_ChemComp[i]->second->get_name() ); 

	// tell me how many eta and how many xi we have
	size_t n_eta, n_xi_mob, n_xi_immob, n_eta_mob, n_eta_immob; 
	// get n_eta and n_xi
	n_eta       = this->_ReductionKin->get_n_eta();
	n_eta_mob   = this->_ReductionKin->get_n_eta_mob(); 
	n_eta_immob = n_eta - n_eta_mob; 
	n_xi_mob    = this->_ReductionKin->get_n_xi_mob(); 
	n_xi_immob  = this->_ReductionKin->get_n_xi_immob(); 

	// creating local memory space to store IC and BC
	// initialize eta_mob, 
	for ( i=0; i < n_eta_mob ; i++ )
	{
		MyNodalFunctionScalar* eta_i = new MyNodalFunctionScalar(); 
		eta_i->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0 );
	    _eta_mob.push_back(eta_i); 
	}
	// initialize eta_immob, 
	for ( i=0; i < n_eta_immob; i++ )
	{
		MyNodalFunctionScalar* eta_i = new MyNodalFunctionScalar(); 
	    eta_i->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0 );
		_eta_immob.push_back(eta_i); 
	}
	// initialize xi_mob
	for ( i=0; i < n_xi_mob ; i++ )
	{
		MyNodalFunctionScalar* xi_mob_tmp       = new MyNodalFunctionScalar();  // xi_mob
        MyNodalFunctionScalar* xi_mob_rates_tmp = new MyNodalFunctionScalar();  // xi_mob_rates
		xi_mob_tmp->initialize(       *dis, FemLib::PolynomialOrder::Linear, 0.0  );
        xi_mob_rates_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_mob.push_back(xi_mob_tmp); 
        _xi_mob_rates.push_back(xi_mob_rates_tmp); 
	}
	// initialize drates_dxi
	for ( i=0; i < n_xi_mob ; i++ )
	{
		MyNodalFunctionScalar* drate_dxi_tmp = new MyNodalFunctionScalar();  // drate_dxi instances
		drate_dxi_tmp->initialize(       *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_mob_drates_dxi.push_back(drate_dxi_tmp); 
	}
	// initialize xi_mob
	for ( i=0; i < n_xi_immob ; i++ )
	{
		MyNodalFunctionScalar* xi_immob_tmp       = new MyNodalFunctionScalar();  // xi_immob
        MyNodalFunctionScalar* xi_immob_rates_tmp = new MyNodalFunctionScalar();  // xi_mob_rates
		MyNodalFunctionScalar* xi_immob_new_tmp   = new MyNodalFunctionScalar();  // xi_mob_rates
		xi_immob_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
        xi_immob_rates_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		xi_immob_new_tmp->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_immob.push_back(xi_immob_tmp); 
        _xi_immob_rates.push_back(xi_immob_rates_tmp); 
		_xi_immob_new.push_back(xi_immob_new_tmp); 
	}

	// linear assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyLinearResidualAssemblerType* linear_r_assembler = new MyLinearResidualAssemblerType(_feObjects);
    MyLinearJacobianAssemblerType* linear_j_eqs = new MyLinearJacobianAssemblerType(_feObjects);

	MyNonLinearAssemblerType* non_linear_assembler = new MyNonLinearAssemblerType(_feObjects, this->_ReductionKin );
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this->_ReductionKin );
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this->_ReductionKin, this);

	_local_ode_xi_immob = new Local_ODE_Xi_immob( this->_ReductionKin ); 

	// for the linear transport problem, variables are eta_mobile
	for ( i=0; i < n_eta_mob ; i++ )
	{
		// set up problem
		MyLinearTransportProblemType* linear_problem = new MyLinearTransportProblemType(dis);
		MyLinearEquationType* linear_eqs = linear_problem->createEquation();
		linear_eqs->initialize(linear_assembler, linear_r_assembler, linear_j_eqs);
		linear_problem->setTimeSteppingFunction(*tim);
		
		// set up variables
		// in this case, the variables includes: eta_0, eta_1, eta_2......, 
		// which are the concentrations of all eta_mobile 
		std::stringstream str_tmp;
		str_tmp << "eta_" << i ; 
		linear_problem->addVariable( str_tmp.str() );
		_linear_problems.push_back(linear_problem); 
	}
	
	// define nonlinear problem
	_non_linear_problem = new MyNonLinearReactiveTransportProblemType(dis);  
	_non_linear_eqs     = _non_linear_problem->createEquation(); 
	_non_linear_eqs->initialize( non_linear_assembler, non_linear_r_assembler, non_linear_j_assembler ); 
	_non_linear_problem->setTimeSteppingFunction(*tim); 
	// for nonlinear coupled transport problem, variables are xi_mobile species
	for ( i=0; i < n_xi_mob ; i++ )
	{
		std::stringstream str_tmp;
		str_tmp << "xi_mob_" << i ;
		_non_linear_problem->addVariable( str_tmp.str() );  	
	}

	// reduction problem
	_problem = new MyKinReductionProblemType( dis, _ReductionKin ); 
	_problem->setTimeSteppingFunction(*tim);  // applying the same time stepping function for all linear, non-linear and reduction problems
	// add variables to the KinReduction problem class
	// for the KinReduction problem, variables are the concentrations of all chemical components
	// add all concentrations to discretized memory space
	for ( i=0; i<n_Comp; i++ )
	{
		_problem->addVariable( femData->map_ChemComp[i]->second->get_name() );
	}

	// set IC for concentrations
	NumLib::TXFunctionBuilder f_builder;
	const BaseLib::Options* opICList = option.getSubGroup("ICList");
    for (const BaseLib::Options* opIC = opICList->getFirstSubGroup("IC"); opIC!=0; opIC = opICList->getNextSubGroup())
    {		
		SolutionLib::FemIC* conc_ic = new SolutionLib::FemIC(msh);
        std::string var_name = opIC->getOption("Variable");
        std::string geo_type = opIC->getOption("GeometryType");
        std::string geo_name = opIC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string dis_name = opIC->getOption("DistributionType");
        double dis_v = opIC->getOptionAsNum<double>("DistributionValue");
        NumLib::ITXFunction* f_ic =  f_builder.create(dis_name, dis_v);
        // set IC
		conc_ic->addDistribution( geo_obj, f_ic ); 
		size_t comp_idx = femData->map_ChemComp.find( var_name )->second->getIndex();  
		_problem->getVariable(comp_idx)->setIC( conc_ic ); 
	}

	// backward flushing global vector of concentrations
	MyNodalFunctionScalar* tmp_conc; 
	for (i=0; i < n_Comp; i++)
	{
		SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
	    tmp_conc = new MyNodalFunctionScalar();
		if ( femIC )
		{
			// FemIC vector is not empty
			tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 0.0);
			femIC->setup(*tmp_conc);
		}
		else
		{
			// FemIC vector is empty
			// initialize the vector with zeros
			tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 0.0);
		}	
		_concentrations.push_back( tmp_conc ); 
	}

	// convert IC _concentrations to eta and xi
	convert_conc_to_eta_xi();

	// set IC for eta_mob
	for ( i=0; i < n_eta_mob; i++ )
	{
		SolutionLib::FemIC* eta_ic = new SolutionLib::FemIC(msh);
		eta_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _eta_mob[i]->getDiscreteData() ) ); 
		_linear_problems[i]->getVariable(0)->setIC( eta_ic ); 
	}
	// set IC for xi_mob
	for ( i=0; i < n_xi_mob; i++ )
	{
		SolutionLib::FemIC* xi_mob_ic = new SolutionLib::FemIC(msh); 
		xi_mob_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _xi_mob[i]->getDiscreteData() ) ); 
		_non_linear_problem->getVariable(i)->setIC( xi_mob_ic ); 
	}
	// set BC for concentrations
	const BaseLib::Options* opBCList = option.getSubGroup("BCList");
    for (const BaseLib::Options* opBC = opBCList->getFirstSubGroup("BC"); opBC!=0; opBC = opBCList->getNextSubGroup())
    {
		std::string var_name = opBC->getOption("Variable");
        std::string geo_type = opBC->getOption("GeometryType");
        std::string geo_name = opBC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string dis_name = opBC->getOption("DistributionType");
        double dis_v = opBC->getOptionAsNum<double>("DistributionValue");
        NumLib::ITXFunction* f_bc =  f_builder.create(dis_name, dis_v);
		size_t comp_idx = femData->map_ChemComp.find( var_name )->second->getIndex();
		_problem->getVariable(comp_idx)->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
    }

	// set ST for concentrations
	const BaseLib::Options* opSTList = option.getSubGroup("STList");
    for (const BaseLib::Options* opST = opSTList->getFirstSubGroup("ST"); opST!=0; opST = opSTList->getNextSubGroup())
    {
		std::string var_name = opST->getOption("Variable");
        std::string geo_type = opST->getOption("GeometryType");
        std::string geo_name = opST->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string st_type = opST->getOption("STType");
        std::string dis_name = opST->getOption("DistributionType");
        double dis_v = opST->getOptionAsNum<double>("DistributionValue");
        if (st_type.compare("NEUMANN")==0) {
            dis_v *= -1; // user set inflow as positive sign but internally negative
        }
        NumLib::ITXFunction* f_st =  f_builder.create(dis_name, dis_v);
        if (f_st!=NULL) {
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(msh, geo_obj, f_st);
            }
            size_t comp_idx = femData->map_ChemComp.find( var_name )->second->getIndex();
			_problem->getVariable(comp_idx)->addNeumannBC(femSt);
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }

    // set up linear solution
	for ( i=0; i < n_eta_mob; i++ )
	{
		MyLinearSolutionType* linear_solution = new MyLinearSolutionType( dis, this->_linear_problems[i] ); 
		MyLinearSolver* linear_solver = linear_solution->getLinearEquationSolver();
		const BaseLib::Options* optNum = option.getSubGroup("Numerics");
		linear_solver->setOption(*optNum);
		this->_linear_solutions.push_back( linear_solution ); 
	}
	
	// set up non-linear solution
	// NRIterationStepInitializer
	MyNRIterationStepInitializer* myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler); 
	MyDiscreteNonlinearSolverFactory* myNSolverFactory = new MyDiscreteNonlinearSolverFactory( myNRIterator ); 
	this->_non_linear_solution = new MyNonLinearSolutionType( dis, this->_non_linear_problem, myNSolverFactory ); 
    this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
    this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver(); 
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);

	// set up solution
    _solution = new MyKinReductionSolution(dis, _problem, this, _linear_problems, _linear_solutions, _non_linear_problem, _non_linear_solution);
	
    // set initial output
    // OutputVariableInfo var(this->getOutputParameterName(Concentrations), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    // femData->outController.setOutput(var.name, var); 

    // initial output parameter
	// TODO
    // this->setOutput(Concentrations, Concentrations->getIC());

    return true;
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	size_t i; 
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for ( i=0; i < _linear_problems.size(); i++ ) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
	// set velocity for nonlinear problem as well
	_non_linear_problem->getEquation()->getLinearAssembler()->setVelocity(vel); 
    _non_linear_problem->getEquation()->getResidualAssembler()->setVelocity(vel); 
	_non_linear_problem->getEquation()->getJacobianAssembler()->setVelocity(vel); 

	_non_linear_problem->getEquation()->getLinearAssembler()->set_xi_mob_rates( &_xi_mob_rates ); 
	_non_linear_problem->getEquation()->getResidualAssembler()->set_xi_mob_rates( &_xi_mob_rates ); 
	_non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_mob_rates( &_xi_mob_rates ); 

	_non_linear_problem->getEquation()->getLinearAssembler()->set_xi_immob_rates( &_xi_mob_rates ); 
	_non_linear_problem->getEquation()->getResidualAssembler()->set_xi_immob_rates( &_xi_mob_rates ); 
	_non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_immob_rates( &_xi_mob_rates ); 
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{ 
	// first use the eta and xi values to calculate comp concentrations
	// TODO
	// then update the output
	// TODO
    // setOutput(Concentrations, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    // update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();

	// convert eta and xi back to concentrations
    convert_eta_xi_to_conc(); 

	// set the new output
	for (size_t i=0; i<_concentrations.size(); i++) {
		OutputVariableInfo var1(this->getOutputParameterName(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
        femData->outController.setOutput(var1.name, var1);
    }
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::convert_conc_to_eta_xi(void)
{
	size_t node_idx, i; 
	size_t n_comp, n_eta_mob, n_eta_immob, n_xi_mob, n_xi_immob; 

	n_comp      = this->_ReductionKin->get_n_Comp();
	n_eta_mob   = this->_ReductionKin->get_n_eta_mob(); 
	n_eta_immob = this->_ReductionKin->get_n_eta() - this->_ReductionKin->get_n_eta_mob() ;
	n_xi_mob    = this->_ReductionKin->get_n_xi_mob(); 
	n_xi_immob  = this->_ReductionKin->get_n_xi_immob(); 
	// only when the reduction scheme is fully initialized
	if ( this->_ReductionKin->IsInitialized() )
	{
		// local vectors
		LocalVector loc_eta_mob;
		LocalVector loc_eta_immob;
		LocalVector loc_xi_mob;
		LocalVector loc_xi_immob;
		LocalVector loc_conc; 
		// allocate the memory for local vectors
		loc_eta_mob     = LocalVector::Zero( this->_ReductionKin->get_n_eta_mob() ); 
		loc_eta_immob   = LocalVector::Zero( this->_ReductionKin->get_n_eta() - this->_ReductionKin->get_n_eta_mob() ); 
		loc_xi_mob      = LocalVector::Zero( this->_ReductionKin->get_n_xi_mob() ); 
		loc_xi_immob    = LocalVector::Zero( this->_ReductionKin->get_n_xi_immob() ); 
		loc_conc        = LocalVector::Zero( this->_ReductionKin->get_n_Comp() );

		// for each nodes, 
		for (node_idx=_concentrations[0]->getDiscreteData()->getRangeBegin(); 
			 node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd(); 
			 node_idx++ )
		{
			for (i=0; i < n_comp; i++)
			{
				// gether all the concentrations 
				loc_conc[i] = _concentrations[i]->getValue(node_idx); 
			}  // end of for i
			
			// pass them to the transform function in the reductionKin class
			// and thet the loc_eta_mob, local_eta_immob and local_xi
			this->_ReductionKin->Conc2EtaXi( loc_conc, loc_eta_mob, loc_eta_immob, loc_xi_mob, loc_xi_immob );
			
			// put the local eta and xi into the global vector
			// fill in eta_mob
			for (i=0; i < n_eta_mob; i++)
				this->_eta_mob[i]->setValue(node_idx, loc_eta_mob[i]); 
			// fill in eta_immob
			for (i=0; i < n_eta_immob; i++)
				this->_eta_immob[i]->setValue(node_idx, loc_eta_immob[i]); 
			// fill in xi_mob
			for (i=0; i < n_xi_mob; i++)
				this->_xi_mob[i]->setValue(node_idx, loc_xi_mob[i]); 
			for (i=0; i < n_xi_immob; i++)
				this->_xi_immob[i]->setValue(node_idx, loc_xi_immob[i]); 
		}  // end of for node_idx
	
	}  // end of if _ReductionKin
}


template <class T1, class T2>
void FunctionConcentrations<T1, T2>::convert_eta_xi_to_conc(void)
{
	size_t node_idx, i; 
	size_t n_comp, n_eta_mob, n_eta_immob, n_xi_mob, n_xi_immob; 

	n_comp      = this->_ReductionKin->get_n_Comp();
	n_eta_mob   = this->_ReductionKin->get_n_eta_mob(); 
	n_eta_immob = this->_ReductionKin->get_n_eta() - this->_ReductionKin->get_n_eta_mob() ;
	n_xi_mob    = this->_ReductionKin->get_n_xi_mob(); 
	n_xi_immob    = this->_ReductionKin->get_n_xi_immob(); 
	// only when the reduction scheme is fully initialized
	if ( this->_ReductionKin->IsInitialized() )
	{
		// local vectors
		LocalVector loc_eta_mob;
		LocalVector loc_eta_immob;
		LocalVector loc_xi_mob;
		LocalVector loc_xi_immob;
		LocalVector loc_conc; 
		// allocate the memory for local vectors
		loc_eta_mob     = LocalVector::Zero( this->_ReductionKin->get_n_eta_mob() ); 
		loc_eta_immob   = LocalVector::Zero( this->_ReductionKin->get_n_eta() - this->_ReductionKin->get_n_eta_mob() ); 
		loc_xi_mob      = LocalVector::Zero( this->_ReductionKin->get_n_xi_mob() ); 
		loc_xi_immob    = LocalVector::Zero( this->_ReductionKin->get_n_xi_immob() ); 
		loc_conc        = LocalVector::Zero( this->_ReductionKin->get_n_Comp() );

		// for each nodes, 
		for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin(); 
			 node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd()  ; 
			 node_idx++ )
		{
			// put the local eta and xi into the global vector
			// fill in eta_mob
			for (i=0; i < n_eta_mob; i++)
				loc_eta_mob[i] = this->_eta_mob[i]->getValue(node_idx); 
			// fill in eta_immob
			for (i=0; i < n_eta_immob; i++)
				loc_eta_immob[i] = this->_eta_immob[i]->getValue(node_idx); 
			// fill in xi
			// this->_xi->setNodalValues( &loc_xi, node_idx*n_xi, n_xi ); 
			for (i=0; i < n_xi_mob; i++)
				loc_xi_mob[i] = this->_xi_mob[i]->getValue(node_idx); 
		    for (i=0; i < n_xi_immob; i++)
				loc_xi_immob[i] = this->_xi_immob[i]->getValue(node_idx); 

			// pass them to the transform function in the reductionKin class
			// and thet the loc_eta_mob, local_eta_immob and local_xi
			this->_ReductionKin->EtaXi2Conc(loc_eta_mob, loc_eta_immob, loc_xi_mob, loc_xi_immob, loc_conc);

			for (i=0; i < n_comp; i++)
			{
				// gether all the concentrations 
				_concentrations[i]->setValue(node_idx, loc_conc[i]); 
			}  // end of for i
		}  // end of for node_idx
	}  // end of if _ReductionKin
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::set_eta_mob_node_values( size_t eta_mob_idx, MyNodalFunctionScalar* new_eta_mob_node_values )
{
	this->_eta_mob[eta_mob_idx] = new_eta_mob_node_values; 
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::set_eta_immob_node_values( size_t eta_immob_idx, MyNodalFunctionScalar* new_eta_immob_node_values )
{
	this->_eta_immob[eta_immob_idx] = new_eta_immob_node_values; 
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::set_xi_mob_node_values( size_t xi_mob_idx, MyNodalFunctionScalar* new_xi_mob_node_values )
{
	this->_xi_mob[xi_mob_idx] = new_xi_mob_node_values; 
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::set_xi_immob_node_values( size_t xi_immob_idx, MyNodalFunctionScalar* new_xi_immob_node_values ) 
{
	this->_xi_immob[xi_immob_idx] = new_xi_immob_node_values; 
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::update_node_kin_reaction_rates(void)
{
	size_t node_idx, i; 
	size_t n_comp, n_eta_mob, n_eta_immob, n_xi_mob, n_xi_immob; 

	n_comp      = this->_ReductionKin->get_n_Comp();
	n_eta_mob   = this->_ReductionKin->get_n_eta_mob(); 
	n_eta_immob = this->_ReductionKin->get_n_eta() - this->_ReductionKin->get_n_eta_mob() ;
	n_xi_mob    = this->_ReductionKin->get_n_xi_mob(); 
	n_xi_immob    = this->_ReductionKin->get_n_xi_immob(); 
	// only when the reduction scheme is fully initialized
	if ( this->_ReductionKin->IsInitialized() )
	{
		// local vectors
		LocalVector loc_eta_mob;
		LocalVector loc_eta_immob;
		LocalVector loc_xi_mob;
		LocalVector loc_xi_immob;
		LocalVector loc_conc; 
        LocalVector loc_xi_mob_rates; 
        LocalVector loc_xi_immob_rates; 

		// allocate the memory for local vectors
		loc_eta_mob        = LocalVector::Zero( n_eta_mob ); 
		loc_eta_immob      = LocalVector::Zero( n_eta_immob ); 
		loc_xi_mob         = LocalVector::Zero( n_xi_mob ); 
		loc_xi_immob       = LocalVector::Zero( n_xi_immob ); 
		loc_conc           = LocalVector::Zero( n_comp );
        loc_xi_mob_rates   = LocalVector::Zero( n_xi_mob );
        loc_xi_immob_rates = LocalVector::Zero( n_xi_immob );

		// for each nodes, 
		for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin(); 
			 node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd()  ; 
			 node_idx++ )
		{
			// put the local eta and xi into the global vector
			// fill in eta_mob
			for (i=0; i < n_eta_mob; i++)
				loc_eta_mob[i] = this->_eta_mob[i]->getValue(node_idx); 
			// fill in eta_immob
			for (i=0; i < n_eta_immob; i++)
				loc_eta_immob[i] = this->_eta_immob[i]->getValue(node_idx); 
			// fill in xi
			// this->_xi->setNodalValues( &loc_xi, node_idx*n_xi, n_xi ); 
			for (i=0; i < n_xi_mob; i++)
				loc_xi_mob[i] = this->_xi_mob[i]->getValue(node_idx); 
		    for (i=0; i < n_xi_immob; i++)
				loc_xi_immob[i] = this->_xi_immob[i]->getValue(node_idx); 

            // calculate rates; 
            this->_ReductionKin->Calc_Xi_Rate( loc_eta_mob, 
                                               loc_eta_immob, 
                                               loc_xi_mob, 
                                               loc_xi_immob, 
                                               loc_xi_mob_rates,
                                               loc_xi_immob_rates );

            // write the rates into the global nodal vector
            for (i=0; i < n_xi_mob; i++)
                this->_xi_mob_rates[i]->setValue( node_idx, loc_xi_mob_rates[i] ); 

            for (i=0; i < n_xi_immob; i++)
                this->_xi_mob_rates[i]->setValue( node_idx, loc_xi_immob_rates[i] );

		}  // end of for node_idx
	}  // end of if _ReductionKin
}


template <class T1, class T2>
void FunctionConcentrations<T1, T2>::update_node_kin_reaction_drates_dxi(void)
{
	const double epsilon = 1.0e-6;

    size_t node_idx, i, j;
	size_t n_eta_mob, n_eta_immob, n_xi_mob, n_xi_immob; 

	LocalVector loc_eta_mob;
	LocalVector loc_eta_immob;
	LocalVector loc_xi_mob;
	LocalVector loc_xi_mob_tmp;  // used for xi increment
	LocalVector loc_xi_immob;
	LocalVector loc_xi_mob_rates_base; 
    LocalVector loc_xi_mob_rates_new; 
    LocalVector loc_xi_immob_rates; 

	n_eta_mob   = this->_ReductionKin->get_n_eta_mob(); 
	n_eta_immob = this->_ReductionKin->get_n_eta_immob();
	n_xi_mob    = this->_ReductionKin->get_n_xi_mob(); 
	n_xi_immob    = this->_ReductionKin->get_n_xi_immob(); 

	// initialize the local vector
	loc_eta_mob           = LocalVector::Zero( n_eta_mob ); 
	loc_eta_immob         = LocalVector::Zero( n_eta_immob ); 
	loc_xi_mob            = LocalVector::Zero( n_xi_mob ); 
	loc_xi_immob          = LocalVector::Zero( n_xi_immob );
    loc_xi_mob_rates_base = LocalVector::Zero( n_xi_mob );
    loc_xi_mob_rates_new  = LocalVector::Zero( n_xi_mob );
    loc_xi_immob_rates    = LocalVector::Zero( n_xi_immob );

	// loop over all nodes
	for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin(); 
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd()  ; 
		 node_idx++)
	{
        // read the local values
        for (i=0; i < n_eta_mob; i++)
			loc_eta_mob[i] = this->_eta_mob[i]->getValue(node_idx); 
		// fill in eta_immob
		for (i=0; i < n_eta_immob; i++)
			loc_eta_immob[i] = this->_eta_immob[i]->getValue(node_idx); 
		for (i=0; i < n_xi_mob; i++)
			loc_xi_mob[i] = this->_xi_mob[i]->getValue(node_idx); 
	    for (i=0; i < n_xi_immob; i++)
			loc_xi_immob[i] = this->_xi_immob[i]->getValue(node_idx); 
        
        // calculate rates; 
        this->_ReductionKin->Calc_Xi_Rate( loc_eta_mob, 
                                           loc_eta_immob, 
                                           loc_xi_mob, 
                                           loc_xi_immob, 
                                           loc_xi_mob_rates_base,
                                           loc_xi_immob_rates );
		// loop over all xi
		for (i=0; i < n_xi_mob; i++)
		{
            // loop over each xi_mob, 
            for (j=0; j < n_xi_mob; j++)
            {
                // get a clean copy of the origiinal xi_mob
                loc_xi_mob_tmp = loc_xi_mob; 

    			// give an increment to the xi value
                loc_xi_mob_tmp(j) += epsilon * loc_xi_mob_tmp(j);

    			// calculate the new rate
                this->_ReductionKin->Calc_Xi_Rate( loc_eta_mob, 
                                                   loc_eta_immob, 
                                                   loc_xi_mob, 
                                                   loc_xi_immob, 
                                                   loc_xi_mob_rates_new,
                                                   loc_xi_immob_rates );

    			// divide the rate value by delta_xi to get derivative
                _xi_mob_drates_dxi[i*n_xi_mob+j]->setValue( node_idx, ( loc_xi_mob_rates_new(i) - loc_xi_mob_rates_base(i) ) / epsilon * loc_xi_mob_tmp(j) ); 
	
            }  // end of for j
		}  // end of for i
	}  // end of for node_idx

    // setting the rates 
    _non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_mob_drate_dxi( &_xi_mob_drates_dxi );
}


template <class T1, class T2>
void FunctionConcentrations<T1, T2>::calc_nodal_xi_immob_ode(double dt)
{
	size_t i; 
	
	// initialize the local vector
	MathLib::LocalVector loc_eta_mob; 
	MathLib::LocalVector loc_eta_immob; 
	MathLib::LocalVector loc_xi_mob; 
	MathLib::LocalVector loc_xi_immob;
	MathLib::LocalVector loc_xi_immob_new; 

	// parameter initialization
	size_t node_idx; 
	size_t n_eta_mob   = this->_ReductionKin->get_n_eta_mob(); 
	size_t n_eta_immob = this->_ReductionKin->get_n_eta_immob();
	size_t n_xi_mob    = this->_ReductionKin->get_n_xi_mob(); 
	size_t n_xi_immob  = this->_ReductionKin->get_n_xi_immob(); 

	// initialize the local vector
	loc_eta_mob           = LocalVector::Zero( n_eta_mob ); 
	loc_eta_immob         = LocalVector::Zero( n_eta_immob ); 
	loc_xi_mob            = LocalVector::Zero( n_xi_mob ); 
	loc_xi_immob          = LocalVector::Zero( n_xi_immob );
	loc_xi_immob_new      = LocalVector::Zero( n_xi_immob );

	// signal, solving local ODEs of xi_immob
	INFO("--Solving local ODE problem of xi_immob...");

	// initialize the ODE Runge-Kutta solution class
	MathLib::RungeKutta4<Local_ODE_Xi_immob, MathLib::LocalVector>* rk4 = new MathLib::RungeKutta4<Local_ODE_Xi_immob, MathLib::LocalVector>(); 

	// loop over all the nodes
	for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin(); 
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd(); 
		 node_idx++ )
	{
		// skip the boundary nodes
		if ( ! this->_solution->isBCNode(node_idx) )
		{
			// on each node, get the right start value
			// get the right set of eta and xi
			for (i=0; i < n_eta_mob; i++)
				loc_eta_mob[i] = this->_eta_mob[i]->getValue(node_idx); 
			// fill in eta_immob
			for (i=0; i < n_eta_immob; i++)
				loc_eta_immob[i] = this->_eta_immob[i]->getValue(node_idx); 
			for (i=0; i < n_xi_mob; i++)
				loc_xi_mob[i] = this->_xi_mob[i]->getValue(node_idx); 
			for (i=0; i < n_xi_immob; i++)
				loc_xi_immob[i] = this->_xi_immob[i]->getValue(node_idx); 
		
			// get the right reference values to ODE RHS function
			this->_local_ode_xi_immob->update_eta_xi( loc_eta_mob, loc_eta_immob, loc_xi_mob, loc_xi_immob); 

			// solve the local ODE problem using RK
			rk4->solve( *_local_ode_xi_immob, 0.0, dt, loc_xi_immob, loc_xi_immob_new); 
		
			// collect the new xi_immob_new
			for (i=0; i < n_xi_immob; i++)
				_xi_immob_new[i]->setValue(node_idx, loc_xi_immob_new[i]); 
		} // end of if

	}  // end of for node_idx

}
