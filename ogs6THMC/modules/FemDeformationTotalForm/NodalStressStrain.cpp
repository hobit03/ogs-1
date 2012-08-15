/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodalStressStrain.cpp
 *
 * Created on 2012-08-14 by Norihiro Watanabe
 */

#include "NodalStressStrain.h"

#include "logog.hpp"

#include "FemLib/Post/Extrapolation.h"

#include "Ogs6FemData.h"
#include "FemLinearElasticTools.h"


size_t FunctionNodalStressStrain::getNumberOfStrainComponents() const
{
    MeshLib::IMesh* msh = _dis->getMesh();
    const size_t dim = msh->getDimension();
    return (dim==2 ? 4 : 6);
}

bool FunctionNodalStressStrain::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOption<size_t>("MeshID");
    _dis = femData->list_dis_sys[msh_id];
    const size_t n_strain_components = getNumberOfStrainComponents();

    // create strain, stress vectors
    NumLib::LocalVector v0(n_strain_components);
    v0 *= .0;

    // set initial output
    _nodal_strain = new FemLib::FemNodalFunctionVector(*_dis, FemLib::PolynomialOrder::Linear, v0);
    _nodal_stress = new FemLib::FemNodalFunctionVector(*_dis, FemLib::PolynomialOrder::Linear, v0);
    for (size_t i=0; i<n_strain_components; i++) {
        _vec_nodal_strain_components.push_back(new NodalPointScalarWrapper(_nodal_strain, i));
        _vec_nodal_stress_components.push_back(new NodalPointScalarWrapper(_nodal_stress, i));
    }
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(getStrainComponentName(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(getStressComponentName(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }

    // initial output parameter
    setOutput(NodStrain, _nodal_strain);
    setOutput(NodStress, _nodal_stress);

    return true;
}

void FunctionNodalStressStrain::accept(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t n_strain_components = getNumberOfStrainComponents();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(getStrainComponentName(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(getStressComponentName(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }
};

int FunctionNodalStressStrain::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Solving NODAL_STRESS_STRAIN...");

    FemLib::FEMIntegrationPointFunctionVector* strain = (FemLib::FEMIntegrationPointFunctionVector*)getInput(GpStrain);
    FemLib::FEMIntegrationPointFunctionVector* stress = (FemLib::FEMIntegrationPointFunctionVector*)getInput(GpStress);

    FemLib::FemExtrapolationAverage<NumLib::LocalVector> extrapo;
    extrapo.extrapolate(*strain, *_nodal_strain);
    extrapo.extrapolate(*stress, *_nodal_stress);

    setOutput(NodStrain, _nodal_strain);
    setOutput(NodStress, _nodal_stress);

    return 0;
}
