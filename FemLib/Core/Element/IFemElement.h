
#pragma once

#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IElement.h"

#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Core/ShapeFunction/ShapeFunction.h"
#include "FiniteElementType.h"

namespace FemLib
{

typedef MathLib::Matrix<double> LocalMatrix;

/**
 * \brief IFiniteElement class is an interface to all kinds of finite element classes. 
 *
 * Roles of the class is
 * - assembly of typical local matrices, e.g. N^T N, dN^T dN
 * - interpolation, extrapolation
 * - integration over domain, boundary
 */
class IFiniteElement
{
public:
	virtual ~IFiniteElement() {};

    /// setup object for given mesh elements
    virtual void configure(MeshLib::IElement &e ) = 0;

    /// return finite element type
    virtual FiniteElementType::type getFeType() const = 0;

    /// return this mesh element
    virtual MeshLib::IElement* getElement() const = 0;

    /// return the number of variables
    virtual size_t getNumberOfVariables() const = 0;

    /// compute basis functions \f$ \mathbf{N}_e \f$, \f$ {\nabla}_x \mathbf{N}_e \f$ at the given point
    virtual void computeBasisFunctions(const double *x_ref) = 0;

    /// compute real coordinates from the given position in reference coordinates
    virtual void getRealCoordinates(double* x_real) = 0;

    /// get evaluated basis functions \f$ \mathbf{N}_e \f$. 
    virtual LocalMatrix* getBasisFunction() = 0;

    /// get evaluated gradient of basis functions \f$ {\nabla}_x \mathbf{N}_e \f$. 
    virtual LocalMatrix* getGradBasisFunction() = 0;

    /// make interpolation from nodal values \f$ u^h(\mathbf{x}) = \mathbf{N}_e (\mathbf x) \mathbf{u}_e \f$
    virtual double interpolate(double *pt, double *nodal_values) = 0;
//    /// compute an matrix \f$ \mathbf{M}_e = \int_{\Omega_e} {\mathbf{N}_e^*}^T f(\mathbf x) \mathbf{N}_e d\Omega \f$
//    virtual void integrateWxN(MathLib::SpatialFunctionScalar*, LocalMatrix &) = 0;
//    /// compute an matrix \f$ \mathbf{M}_e = \int_{\Omega_e} {\mathbf{N}_e^*}^T \mathbf{f}(\mathbf x) \nabla \mathbf{N}_e d\Omega \f$
//    virtual void integrateWxDN(MathLib::SpatialFunctionVector*, LocalMatrix &) = 0;
//    /// compute an matrix \f$ \mathbf{M}_e = \int_{\Omega_e} {\nabla \mathbf{N}_e^*}^T f(\mathbf x) \nabla \mathbf{N}_e d\Omega \f$
//    virtual void integrateDWxDN( MathLib::SpatialFunctionScalar* f, LocalMatrix &) = 0;

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(size_t igp, double f, LocalMatrix &mat) = 0;

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(size_t igp, MathLib::Vector &f, LocalMatrix &mat) = 0;

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(size_t igp, double f, LocalMatrix &mat) = 0;

    /// get the integration method
    virtual IFemNumericalIntegration* getIntegrationMethod() const = 0;

    virtual void extrapolate(const std::vector<MathLib::Vector> &gp_values, std::vector<MathLib::Vector> &nodal_values) = 0;
};


} //end

