/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file XFEM_EXAMPLE_CRACK1.cpp
 *
 * Created on 2012-09-20 by Norihiro Watanabe
 */

#include "XFEM_EXAMPLE_CRACK1.h"

#include <cmath>
#include <limits>
#include <vector>
#include <set>

#include "logog.hpp"

#include "MathLib/LinAlg/LinearEquation/EigenDenseLinearEquation.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
//#include "FemLib/Function/FemNodalFunction.h"
//#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
//#include "SolutionLib/Fem/FemSourceTerm.h"
//#include "MaterialLib/PorousMedia.h"
//#include "MaterialLib/Solid.h"
#include "Ogs6FemData.h"

#include "../FemDeformationTotalForm/FemLinearElasticTools.h"
#include "matlab_functions.h"

namespace xfem
{

//template <class T1, class T2>
bool FunctionXFEM_EXAMPLE_CRACK1::initialize(const BaseLib::Options &option)
{
    option.printout(std::cout);

    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
//    _msh = MeshLib::MeshGenerator::generateRegularQuadMesh(2., 3, -1., -1., .0);
//    INFO("Mesh:");
//    INFO("* Nr. of nodes = %d", _msh->getNumberOfNodes());
//    INFO("* Nr. of elements = %d", _msh->getNumberOfElements());
    _msh = femData->list_mesh[msh_id];
    _dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(_msh);
    _feObjects = new FemLib::LagrangianFeObjectContainer(*_msh);

    _displacement = new MyNodalFunctionVector();
    NumLib::LocalVector tmp_u0 = NumLib::LocalVector::Zero(2);
    _displacement->initialize(*_dis, FemLib::PolynomialOrder::Linear, tmp_u0);
    _exact_displacement = new MyNodalFunctionVector();
    _exact_displacement->initialize(*_dis, FemLib::PolynomialOrder::Linear, tmp_u0);

    // set initial output
    OutputVariableInfo var("DISPLACEMENT", OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    OutputVariableInfo var_eu("EXACT_U", OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _exact_displacement);
    femData->outController.setOutput(var_eu.name, var_eu);
//    for (size_t i=0; i<_vec_u_components.size(); i++) {
//        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }

//    // initial output parameter
//    this->setOutput(Displacement, _displacement);


    return true;
}


//template <class T1, class T2>
int FunctionXFEM_EXAMPLE_CRACK1::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("-> solve %s", this->getProcessName().c_str());
    const size_t NodeNum = _msh->getNumberOfNodes();
    const size_t ElemNum = _msh->getNumberOfElements();

    // define test case parameters
    const double k1 = 1.0;
    const double EE = 10000.;
    const double nu = 0.3;
    const double fx = 0, fy = 0;
    const double mu = EE/(2.*(1.+nu));
    //lambda = EE*nu/((1+nu)*(1-2*nu)); % plane strain
    //kappa = 3-4*nu; % plane strain
    const double lambda = EE*nu/(1.-nu*nu); // plane stress
    const double kappa  = (3.-nu)/(1.+nu); // plane stress

    // get level-set function
    std::vector<double> ff(NodeNum);
    for (size_t i=0; i<NodeNum; i++)
        ff[i] = _msh->getNodeCoordinatesRef(i)->getData()[1];

    // Get exact solution at the nodes.
    std::vector<double> uuExact(NodeNum, .0), vvExact(NodeNum, .0);
    for (size_t i=0; i<NodeNum; i++) {
        const GeoLib::Point *pt = _msh->getNodeCoordinatesRef(i);
        exactSol_Mode1((*pt)[0], (*pt)[1], k1, kappa, mu, lambda, uuExact[i], vvExact[i]);
        _exact_displacement->getValue(i)(0) = uuExact[i];
        _exact_displacement->getValue(i)(1) = vvExact[i];
    }

    // define Dirichlet BC
    std::vector<size_t> Bound;
    searchMinMaxNodes(*_msh, Bound);
    std::vector<size_t> uDirNodes(Bound), vDirNodes(Bound);
    for (size_t i=0; i<vDirNodes.size(); i++)
        vDirNodes[i] += NodeNum;
    std::vector<double> uDirValues(Bound.size()), vDirValues(Bound.size());
    for (size_t i=0; i<Bound.size(); i++) {
        uDirValues[i] = uuExact[Bound[i]];
        vDirValues[i] = vvExact[Bound[i]];
    }

    // Get enriched elements and nodes.
    std::vector<size_t> ElemsEnriched, NodesEnriched;
    getEnrichedNodesElems(*_msh, ff, ElemsEnriched, NodesEnriched);
    std::set<size_t> SetNodesEnriched;
    SetNodesEnriched.insert(NodesEnriched.begin(), NodesEnriched.end());

    // initialize LinearEQS
    MathLib::EigenDenseLinearEquation leqs;
    leqs.create(4*NodeNum);
    INFO("* Nr. of DoFs = %d", 4*NodeNum);

    // domain integration
    INFO("* start Domain integration");
    for (size_t i=0; i<ElemNum; i++) {

        MeshLib::IElement* e = _msh->getElement(i);
        const size_t n_ele_nodes = e->getNumberOfNodes();
        NumLib::LocalVector Nodes(n_ele_nodes);
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(*e);
        NumLib::LocalVector xxElem(n_ele_nodes), yyElem(n_ele_nodes);
        NumLib::LocalVector ffEle(n_ele_nodes);
        for (size_t j=0; j<n_ele_nodes; j++) {
            Nodes(j) = e->getNodeID(j);
            xxElem(j) = _msh->getNodeCoordinatesRef(e->getNodeID(j))->getData()[0];
            yyElem(j) = _msh->getNodeCoordinatesRef(e->getNodeID(j))->getData()[1];
            ffEle(j) = ff[e->getNodeID(j)];
        }

        // activate nodes are enriched
        NumLib::LocalVector NodesAct(n_ele_nodes);
        for (size_t i=0; i<n_ele_nodes; i++) {
            NodesAct(i) = (SetNodesEnriched.count(Nodes(i))>0 ? 1 : 0);
        }

        // set integration points in the reference element
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        const size_t nQnQ = q->getNumberOfSamplingPoints();
        std::vector<GeoLib::Point> vec_int_ref_xx(nQnQ);
        std::vector<double> vec_int_ref_w(nQnQ);
        for (size_t j=0; j<nQnQ; j++) {
            q->getSamplingPoint(j, vec_int_ref_xx[j].getData());
            vec_int_ref_w[j] = q->getWeight(j);
        }
        NumLib::LocalVector xxIntRef, yyIntRef, wwIntRef;
        IntPoints2DLevelSet(ffEle, vec_int_ref_xx, vec_int_ref_w, nQnQ, xxIntRef, yyIntRef, wwIntRef);
        const size_t Curr_nQ = xxIntRef.rows();

        // get shape functions
        NumLib::LocalMatrix N, dNdx, dNdy;
        NumLib::LocalMatrix M, dMdx, dMdy;
        NumLib::LocalVector xxInt, yyInt, wwInt;
        NumLib::LocalVector ffInt;
        ShapeFctsXFEMSign(
                xxElem, yyElem, ffEle, NodesAct, xxIntRef, yyIntRef, wwIntRef,
                Curr_nQ,
                N, dNdx, dNdy, M, dMdx, dMdy, xxInt, yyInt, wwInt, ffInt);

        // integrate
        BuildMatRhs_Hooke(
                N, dNdx, dNdy, M, dMdx, dMdy,
                xxInt, yyInt, wwInt, ffInt, Nodes,
                lambda, lambda, mu, mu, fx, fy, Curr_nQ, NodeNum,
                leqs);
    }

    // Insert Dirichlet BCs.
    INFO("* insert Dirichlet BCs.");
    leqs.setKnownX(uDirNodes, uDirValues);
    leqs.setKnownX(vDirNodes, vDirValues);
    std::vector<size_t> uNonEnrichedNodes, vNonEnrichedNodes;
    for (size_t i=0; i<NodeNum; i++) {
        if (SetNodesEnriched.count(i) ==0) {
            uNonEnrichedNodes.push_back(i+NodeNum*2);
            vNonEnrichedNodes.push_back(i+NodeNum*3);
        }
    }
    std::vector<double> zeroEnrichedValue(uNonEnrichedNodes.size(), .0);
    leqs.setKnownX(uNonEnrichedNodes, zeroEnrichedValue);
    leqs.setKnownX(vNonEnrichedNodes, zeroEnrichedValue);


//    // Reduce system of equations.
//    Pos = [[1:1:2*NodeNum]'; NodesEnriched+2*NodeNum; NodesEnriched+3*NodeNum];
//    MAT = MAT(Pos, Pos);
//    RHS = RHS(Pos);

    //disp(sprintf('Condition number of final system         : %15.5e', condest(MAT)))

    // Solve system of equations for solution.
    INFO("* solve system of equations");
    leqs.solve();

    double *x = leqs.getX();
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = x[i];
        _displacement->getValue(i)(1) = x[i+NodeNum];
    }


//    std::cout << "x = ";
//    for (size_t i=0; i<leqs.getDimension(); i++)
//        std::cout << x[i] << " ";
//    std::cout << std::endl;
//
//    leqs.printout(std::cout);
//    uuTotal = Sol(1:NodeNum);
//    vvTotal = Sol(NodeNum+1:2*NodeNum);

    return 0;
}

//template <class T1, class T2>
void FunctionXFEM_EXAMPLE_CRACK1::accept(const NumLib::TimeStep &/*time*/)
{
//    //update data for output
//    const size_t n_strain_components = getNumberOfStrainComponents();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
//    for (size_t i=0; i<n_strain_components; i++) {
//        OutputVariableInfo var1(this->getOutputParameterName(NodStrain) + getStressStrainComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_strain_components[i]);
//        femData->outController.setOutput(var1.name, var1);
//        OutputVariableInfo var2(this->getOutputParameterName(NodStress) + getStressStrainComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_stress_components[i]);
//        femData->outController.setOutput(var2.name, var2);
//    }
    OutputVariableInfo var("DISPLACEMENT", OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    OutputVariableInfo var_eu("EXACT_U", OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _exact_displacement);
    femData->outController.setOutput(var_eu.name, var_eu);
}

//template <class T1, class T2>
void FunctionXFEM_EXAMPLE_CRACK1::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
//    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
//        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
//        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
//    }
    //setOutput(Displacement, _displacement);

    //calculateStressStrain();
}

//template <class T1, class T2>
void FunctionXFEM_EXAMPLE_CRACK1::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
};

}
