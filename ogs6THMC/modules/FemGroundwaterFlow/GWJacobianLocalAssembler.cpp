
#include "GWJacobianLocalAssembler.h"

#include "NumLib/Function/TXFunction.h"
#include "Ogs6FemData.h"

void GroundwaterFlowJacobianLocalAssembler::assembly(const NumLib::TimeStep &ts, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ)
{
    FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
    size_t mat_id = e.getGroupID();
    MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
    NumLib::LocalMatrix localM(localJ.rows(), localJ.cols());
    NumLib::LocalMatrix localK(localJ.rows(), localJ.cols());
    localM *= .0;
    localK *= .0;

    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);

        NumLib::LocalMatrix k;
        pm->hydraulic_conductivity->eval(real_x, k);
        NumLib::LocalMatrix s;
        pm->storage->eval(real_x, s);

        fe->integrateWxN(j, s, localM);
        fe->integrateDWxDN(j, k, localK);
    }

    double euler_theta = 1.0; //TODO where to get this parameter?
    localJ.noalias() = 1./ts.getTimeStepSize() * localM + euler_theta *localK;
}

