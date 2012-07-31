
#pragma once

#include <vector>

#include "MeshLib/Core/IElement.h"
#include "NumLib/Function/Function.h"

#include "FemLib/Core/DataType.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Core/ShapeFunction/ShapeFunction.h"
#include "FiniteElementType.h"

namespace FemLib
{

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

    /// get evaluated determinant of the Jacobian matrix for coordinate transformation
    virtual double getDetJ() const = 0;

    /// make interpolation from nodal values \f$ u^h(\mathbf{x}) = \mathbf{N}_e (\mathbf x) \mathbf{u}_e \f$
    virtual double interpolate(double *pt, double *nodal_values) = 0;

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(size_t igp, LocalMatrix &f, LocalMatrix &mat) = 0;

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(size_t igp, LocalMatrix &f, LocalMatrix &mat) = 0;

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(size_t igp, LocalMatrix &f, LocalMatrix &mat) = 0;

    /// get the integration method
    virtual IFemNumericalIntegration* getIntegrationMethod() const = 0;

    virtual void extrapolate(const std::vector<LocalVector> &gp_values, std::vector<LocalVector> &nodal_values) = 0;
};


} //end

