/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticLinearLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "LinearElasticLinearLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"

#include "FemLinearElasticTools.h"
#include "Ogs6FemData.h"

void FemLinearElasticLinearLocalAssembler::assembly(const NumLib::TimeStep &/*time*/,  const MeshLib::IElement &e, const LocalVectorType &/*local_u_n1*/, const LocalVectorType &/*local_u_n*/, NumLib::LocalEquation &eqs)
{
    FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
    size_t mat_id = e.getGroupID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    //MaterialLib::PorousMedia* pm = femData->list_pm[mat_id];
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];

    const size_t dim = e.getDimension();
    const size_t n_strain_components = (dim==2 ? 4 : 6);
    const size_t nnodes = e.getNumberOfNodes();
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // set D
    LocalMatrixType matD = LocalMatrixType::Zero(n_strain_components, n_strain_components);
    //matD *= .0;
    NumLib::LocalMatrix nv(1,1);
    NumLib::LocalMatrix E(1,1);
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);
    double Lambda, G, K;
    MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
    MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);

    // body force
    LocalVectorType body_force(dim);
    body_force *= .0;
    bool hasGravity = false;
    if (hasGravity) {
        double rho_s = .0;
        solidphase->density->eval(e_pos, rho_s);
        body_force[dim-1] = rho_s * 9.81;
    }

    //
    LocalMatrixType &localK = *eqs.getA();
    LocalVectorType &localRHS = *eqs.getRHSAsVec();
    LocalMatrixType matB(n_strain_components, nnodes*dim);
    LocalMatrixType matN(dim, nnodes*dim);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);
        double fac = fe->getDetJ() * q->getWeight(j);


        // set N,B
        LocalMatrixType &N = *fe->getBasisFunction();
        LocalMatrixType &dN = *fe->getGradBasisFunction();
        setNu_Matrix(dim, nnodes, N, matN);
        setB_Matrix(dim, nnodes, dN, matB);

        // K += B^T * D * B
        localK.noalias() += fac * matB.transpose() * matD * matB;

        // RHS
        localRHS.noalias() += fac * matN.transpose() * body_force;
    }
}

