
#include "FemExtrapolation.h"
#include "IFemElement.h"

namespace FemLib
{
// average
void FeExtrapolationAverage::extrapolate(IFiniteElement &fe, const std::vector<MathLib::Vector2D> &gp_values, std::vector<MathLib::Vector2D> &nodal_values)
{
    extrapolate<MathLib::Vector2D>(fe, gp_values, nodal_values);
}

template<typename Tvalue>
void FeExtrapolationAverage::extrapolate(IFiniteElement &fe, const std::vector<Tvalue> &gp_values, std::vector<Tvalue> &nodal_values)
{
    Tvalue ele_avg;
    ele_avg = .0;
    for (size_t j=0; j<gp_values.size(); j++)
        ele_avg += gp_values[j];
    ele_avg /= gp_values.size();

    std::fill(nodal_values.begin(), nodal_values.end(), ele_avg);
}

// gauss linear
void FeExtrapolationGaussLinear::extrapolate(IFiniteElement &fe, const std::vector<MathLib::Vector2D> &gp_values, std::vector<MathLib::Vector2D> &nodal_values)
{
    extrapolate<MathLib::Vector2D>(fe, gp_values, nodal_values);
}

template<typename Tvalue>
void FeExtrapolationGaussLinear::extrapolate(IFiniteElement &fe, const std::vector<Tvalue> &gp_values, std::vector<Tvalue> &nodal_values)
{
    const MeshLib::IElement* e = fe.getElement();
    FemIntegrationGaussBase* gauss = static_cast<FemIntegrationGaussBase*>(fe.getIntegrationMethod());
    const size_t e_nnodes = e->getNumberOfNodes();
    // reorder gp values
    std::vector<Tvalue> reordered_gp_values(gp_values.size());
    switch (e->getElementType()) {
        case MeshLib::ElementType::QUAD:
        case MeshLib::ElementType::HEXAHEDRON:
            for (size_t j=0; j<gp_values.size(); j++) {
                size_t nod_id = getLocalIndex(*e, *gauss, j);
                reordered_gp_values[nod_id] = gp_values[j];
            }
            break;
        default:
            reordered_gp_values.assign(gp_values.begin(), gp_values.end());
    }

    // calculate Xi_p
    double Xi_p = calcXi_p(*e, *gauss);

    // 
    size_t i_s = 0;
    size_t i_e = e_nnodes;
    size_t ish = 0;
    if (e->getElementType()==MeshLib::ElementType::TETRAHEDRON) {
        i_s = 1;
        i_e = e_nnodes + 1;
        ish = 1;
    }

    // 
    double x[3] = {};
    for (size_t j=0; j<e_nnodes; j++) {
        getExtropoGaussPoints(*e, j, Xi_p, x);
        fe.computeBasisFunctions(x);
        MathLib::Matrix<double> *N = fe.getBasisFunction();
        Tvalue EV(.0);
        for(size_t k=i_s; k<i_e; k++)
            EV += reordered_gp_values[k] * (*N)(0,k - ish);
        nodal_values[j] = EV;
    }
}

double FeExtrapolationGaussLinear::calcXi_p(const MeshLib::IElement& e, FemIntegrationGaussBase& gauss)
{
    double Xi_p = 0.0;
    if (e.getElementType() == MeshLib::ElementType::QUAD || e.getElementType() == MeshLib::ElementType::HEXAHEDRON) {
        double r = .0;
        const size_t nGauss = gauss.getSamplingLevel();
        for (size_t gp=0; gp<nGauss; gp++) {
            r = MathLib::GaussLegendre::getPoint(nGauss, gp);
            if (fabs(r) > Xi_p)
                Xi_p = fabs(r);
        }
        r = 1.0 / Xi_p;
        Xi_p = r;
    }

    return Xi_p;
}

int FeExtrapolationGaussLinear::getLocalIndex(const MeshLib::IElement& e, FemIntegrationGaussBase& gauss, size_t igp)
{
    int LoIndex = -1;
    double r,s,t;
    size_t nGauss = gauss.getSamplingLevel();
    size_t gp_r = 0;
    size_t gp_s = 0;
    size_t gp_t = 0;
    const double MKleinsteZahl = std::numeric_limits<double>::epsilon();

    switch (e.getElementType()) {
    case MeshLib::ElementType::QUAD:
        gp_r = igp/nGauss;
        gp_s = igp%nGauss;
        gp_t = 0;
        r = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
        s = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
        if (r > 0.0 && s > 0.0)
            LoIndex = 0;
        else if (r < 0.0 && s > 0.0)
            LoIndex = 1;
        else if (r < 0.0 && s < 0.0)
            LoIndex = 2;
        else if (r > 0.0 && s < 0.0)
            LoIndex = 3;
        else if (fabs(r) < MKleinsteZahl && s > 0.0)
            LoIndex = 4;
        else if (r < 0.0 && fabs(s) < MKleinsteZahl)
            LoIndex = 5;
        else if (fabs(r) < MKleinsteZahl && s < 0.0)
            LoIndex = 6;
        else if (r > 0.0 && fabs(s) < MKleinsteZahl)
            LoIndex = 7;
        else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
            LoIndex = 8;
        break;
    case MeshLib::ElementType::HEXAHEDRON:
        gp_r = igp/(nGauss * nGauss);
        gp_s = igp%(nGauss * nGauss);
        gp_t = gp_s%nGauss;
        r = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
        s = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
        t = MathLib::GaussLegendre::getPoint(nGauss, gp_t);

        if (t > 0.0) {
            if (r > 0.0 && s > 0.0)
                LoIndex = 0;
            else if (r < 0.0 && s > 0.0)
                LoIndex = 1;
            else if (r < 0.0 && s < 0.0)
                LoIndex = 2;
            else if (r > 0.0 && s < 0.0)
                LoIndex = 3;
            else if (fabs(r) < MKleinsteZahl && s > 0.0)
                LoIndex = 8;
            else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 9;
            else if (fabs(r) < MKleinsteZahl && s < 0.0)
                LoIndex = 10;
            else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 11;
            else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                return -1;
        } else if (fabs(t) < MKleinsteZahl) {
            if (fabs(r) < MKleinsteZahl || fabs(s) < MKleinsteZahl)
                return -1;
            if (r > 0.0 && s > 0.0)
                LoIndex = 16;
            else if (r < 0.0 && s > 0.0)
                LoIndex = 17;
            else if (r < 0.0 && s < 0.0)
                LoIndex = 18;
            else if (r > 0.0 && s < 0.0)
                LoIndex = 19;
        } if (t < 0.0) {
            if (r > 0.0 && s > 0.0)
                LoIndex = 4;
            else if (r < 0.0 && s > 0.0)
                LoIndex = 5;
            else if (r < 0.0 && s < 0.0)
                LoIndex = 6;
            else if (r > 0.0 && s < 0.0)
                LoIndex = 7;
            else if (fabs(r) < MKleinsteZahl && s > 0.0)
                LoIndex = 12;
            else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 13;
            else if (fabs(r) < MKleinsteZahl && s < 0.0)
                LoIndex = 14;
            else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 15;
            else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                return -1;
        }
        break;
    default:
        std::cerr << "CElement::GetLocalIndex invalid mesh element type given"
            << std::endl;
    }
    return LoIndex;
}

void FeExtrapolationGaussLinear::getExtropoGaussPoints(const MeshLib::IElement &e, const int i, double Xi_p, double* unit)
{
    int j = 0;
    //
    switch (e.getElementType()) {
    case MeshLib::ElementType::LINE:
        break;
    case MeshLib::ElementType::TRIANGLE:
        switch (i) {
        case 0:
            unit[0] = -0.1666666666667;
            unit[1] = -0.1666666666667;
            break;
        case 1:
            unit[0] = 1.6666666666667;
            unit[1] = -0.1666666666667;
            break;
        case 2:
            unit[0] = -0.1666666666667;
            unit[1] = 1.6666666666667;
            break;
        }
        break;
    case MeshLib::ElementType::QUAD:
        switch (i) {
        case 0:
            unit[0] = Xi_p;
            unit[1] = Xi_p;
            break;
        case 1:
            unit[0] = -Xi_p;
            unit[1] = Xi_p;
            break;
        case 2:
            unit[0] = -Xi_p;
            unit[1] = -Xi_p;
            break;
        case 3:
            unit[0] = Xi_p;
            unit[1] = -Xi_p;
            break;
        }
        break;
    case MeshLib::ElementType::HEXAHEDRON:
        if (i < 4) {
            j = i;
            unit[2] = Xi_p;
        } else  {
            j = i - 4;
            unit[2] = -Xi_p;
        }
        switch (j) {
        case 0:
            unit[0] = Xi_p;
            unit[1] = Xi_p;
            break;
        case 1:
            unit[0] = -Xi_p;
            unit[1] = Xi_p;
            break;
        case 2:
            unit[0] = -Xi_p;
            unit[1] = -Xi_p;
            break;
        case 3:
            unit[0] = Xi_p;
            unit[1] = -Xi_p;
            break;
        }
        break;
    case MeshLib::ElementType::TETRAHEDRON:
        switch (i) {
        case 0:
            unit[0] = -0.166666666666667;
            unit[1] = -0.166666666666667;
            unit[2] = -0.166666666666667;
            break;
        case 1:
            unit[0] = 1.5;
            unit[1] = -0.166666666666667;
            unit[2] = -0.166666666666667;
            break;
        case 2:
            unit[0] = -0.166666666666667;
            unit[1] = 1.5;
            unit[2] = -0.166666666666667;
            break;
        case 3:
            unit[0] = -0.166666666666667;
            unit[1] = -0.166666666666667;
            unit[2] = 1.5;
        }
        break;
    default:
        unit[0] = unit[1] = unit[2] = 0.;
        break;
    }
}


} //end namespace

