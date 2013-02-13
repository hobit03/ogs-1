/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


//#include "Concentration.h"

#include "logog.hpp"

#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"
#include "MathLib/DataType.h"

template <class T1, class T2>
bool FunctionRichards<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    //MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);
	
    // local assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);
	
    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    typename MyProblemType::MyVariable* pressure_w_variable = _problem->addVariable("PRESSURE1");
    FemVariableBuilder varBuilder;
    varBuilder.doit("PRESSURE1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, pressure_w_variable);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MyLinearSolver* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
	
	// set initial output
	CalculateSaturationValuesForOutput(dis);
	GetAndInsertPrimaryVariable(_problem->getDiscreteSystem()->getMesh()->getID());
    GetAndInsertAdditionalOutputVariables(msh_id);    
    return true;
}

template <class T1, class T2>
void FunctionRichards<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
  
}

template <class T1, class T2>
void FunctionRichards<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
	CalculateSaturationValuesForOutput(dis);
    setOutput(Pressure, _solution->getCurrentSolution(0));
	setOutput(Saturation, _saturation_w);
	
}

template <class T1, class T2>
void FunctionRichards<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    //data for output
	GetAndInsertPrimaryVariable(_problem->getDiscreteSystem()->getMesh()->getID());
	GetAndInsertAdditionalOutputVariables(_problem->getDiscreteSystem()->getMesh()->getID());
};


template <class T1, class T2>
void FunctionRichards<T1, T2>::GetAndInsertPrimaryVariable(const size_t msh_id)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var1(getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var1.name, var1); 
};

template <class T1, class T2>
void FunctionRichards<T1, T2>::GetAndInsertAdditionalOutputVariables(const size_t msh_id)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var2(getOutputParameterName(Saturation), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _saturation_w);
    femData->outController.setOutput(var2.name, var2); 
};

template <class T1, class T2>
void FunctionRichards<T1, T2>::CalculateSaturationValuesForOutput(MyDiscreteSystem* dis)
{
	_pressure_w = new MyNodalFunctionScalar();
	_pressure_w = _solution->getCurrentSolution(0);
	_saturation_w = new MyNodalFunctionScalar(); 
	_saturation_w->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);

			// get pointer to corresponding porous media class
			//const size_t n_dim = e.getDimension();
			//size_t mat_id = e.getGroupID();
			//MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
     double Pc = 0.0, Sw=0.0;
	 size_t node_idx; 
	 for (node_idx = _pressure_w->getDiscreteData()->getRangeBegin(); 
			 node_idx < _pressure_w->getDiscreteData()->getRangeEnd(); 
			 node_idx++ )
	 {
	 MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[0];
	 Pc= -1.0 * _pressure_w->getValue(node_idx);
	 Sw= pm->getSwbyPc(Pc);
	 _saturation_w->setValue(node_idx,Sw);	
	 }
};