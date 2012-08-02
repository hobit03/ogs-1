
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"

#include "Ogs6FemData.h"

class MassTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    MassTransportJacobianLocalAssembler(MaterialLib::Compound* cmp, FemLib::LagrangianFeObjectContainer* feObjects)
        : _cmp(cmp), _feObjects(feObjects)
    {
    };

    virtual ~MassTransportJacobianLocalAssembler() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/, NumLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        NumLib::LocalMatrix matM(localJ);
        NumLib::LocalMatrix matDiff(localJ);
        NumLib::LocalMatrix matAdv(localJ);
        NumLib::TXCompositFunction<NumLib::ITXFunction, NumLib::ITXFunction, NumLib::Multiplication> f_diff_poro(*_cmp->molecular_diffusion, *pm->porosity);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            NumLib::LocalMatrix poro(1,1);
            pm->porosity->eval(real_x, poro);
            NumLib::LocalMatrix d_poro(1,1);
            f_diff_poro.eval(real_x, d_poro);
            NumLib::ITXFunction::DataType v;
            _vel->eval(real_x, v);

            fe->integrateWxN(j, poro, matM);
            fe->integrateDWxDN(j, d_poro, matDiff);
            fe->integrateWxDN(j, v, matAdv);
        }

        double dt = time.getTimeStepSize();
        double theta = 1.0;
        matM *= 1.0 / dt;
        matDiff *= theta;
        matAdv *= theta;
        localJ = matM;
        localJ += matDiff;
        localJ += matAdv;

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    MaterialLib::Compound* _cmp;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    NumLib::ITXFunction* _vel;
};
