/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Velocity.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/Vector.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "Tests/Geo/Material/PorousMedia.h"

using namespace NumLib;

namespace Geo
{

class FunctionVelocity
    : public NumLib::AbstractTransientMonolithicSystem
{
    enum In { Head=0 };
    enum Out { Velocity=0 };
public:
    typedef FemLib::FEMIntegrationPointFunctionVector<DiscreteLib::DiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef FemLib::FemNodalFunctionScalar<DiscreteLib::DiscreteSystem>::type MyNodalFunctionScalar;
    
    FunctionVelocity() 
    {
        AbstractTransientMonolithicSystem::resizeInputParameter(1);
        AbstractTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem &dis, PorousMedia &pm)
    {
        _dis = &dis;
        _K = pm.hydraulic_conductivity;
        _vel = new MyIntegrationPointFunctionVector();
        _vel->initialize(&dis);
        //this->setOutput(Velocity, _vel);
    }

    int solveTimeStep(const TimeStep &/*time*/)
    {
        const MeshLib::IMesh *msh = _dis->getMesh();
        MyNodalFunctionScalar *head = (MyNodalFunctionScalar*)getInput(Head);
        MyIntegrationPointFunctionVector *vel = _vel;;

        FemLib::LagrangianFeObjectContainer* feObjects = head->getFeObjectContainer();
        //calculate vel (vel=f(h))
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement* e = msh->getElemenet(i_e);
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            NumLib::LocalVector local_h(e->getNumberOfNodes());
            for (size_t j=0; j<e->getNumberOfNodes(); j++)
                local_h[j] = head->getValue(e->getNodeID(j));
            // for each integration points
            FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
            double r[2] = {};
            const size_t n_gp = integral->getNumberOfSamplingPoints();
            vel->setNumberOfIntegationPoints(i_e, n_gp);
            NumLib::LocalVector xi(e->getNumberOfNodes());
            NumLib::LocalVector yi(e->getNumberOfNodes());
            for (size_t i=0; i<e->getNumberOfNodes(); i++) {
                const GeoLib::Point* pt = msh->getNodeCoordinatesRef(e->getNodeID(i));
                xi[i] = (*pt)[0];
                yi[i] = (*pt)[1];
            }
            for (size_t ip=0; ip<n_gp; ip++) {
                NumLib::LocalVector q(2);
                q[0] = .0;
                q[1] = .0;
                integral->getSamplingPoint(ip, r);
                fe->computeBasisFunctions(r);
                const NumLib::LocalMatrix* dN = fe->getGradBasisFunction();
                NumLib::LocalMatrix* N = fe->getBasisFunction();
                std::vector<double> xx(3, .0);
                NumLib::LocalVector tmp_v;
                tmp_v = (*N) * xi;
                xx[0] = tmp_v[0];
                tmp_v = (*N) * yi;
                xx[1] = tmp_v[0];
                NumLib::TXPosition pos(&xx[0]);

                NumLib::LocalMatrix k;
                _K->eval(pos, k);
                if (k.rows()==1) {
                    q.noalias() = (*dN) * local_h * (-1.0) * k(0,0);
                } else {
                    q.noalias() = (*dN) * k * local_h * (-1.0);
                }
                vel->setIntegrationPointValue(i_e, ip, q);
            }
        }
        setOutput(Velocity, vel);
        return 0;
    }

    double suggestNext(const TimeStep &/*time_current*/)
    {
        return .0;
    }

    bool isAwake(const TimeStep &/*time*/)
    {
        return true;
    }

    void accept(const TimeStep &/*time*/)
    {
        //std::cout << "Velocity=" << std::endl;
        //_vel->printout();
    };

private:
    DiscreteLib::DiscreteSystem* _dis;
    MyIntegrationPointFunctionVector* _vel;
    NumLib::ITXFunction* _K;

    DISALLOW_COPY_AND_ASSIGN(FunctionVelocity);
};

}
