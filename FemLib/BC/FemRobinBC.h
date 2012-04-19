
#pragma once

#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemFunction.h"
#include "IFemBC.h"

namespace FemLib
{

/**
 * \brief Robin (third type) boundary conditions
 *
 * The Robin boundary condition is formulated for the unknown \f$u\f$ as,
 * \f[
 *     \alpha \frac{\partial u}{\partial n} + \beta u = c
 * \f]
 * with non-zero constants \f$\alpha\f$ and \f$\beta\f$ and a given function \f$c\f$.
 * The equation can also be written as
 * \f[
 *     \frac{\partial u}{\partial n} = \frac{c}{\alpha} - \frac{\beta u}{\alpha}.
 * \f]
 * In FEM, this formula can be included into boundary integrations of diffusive terms
 * \f[
 *     \int_{\Gamma} N^* (D \nabla u) \mathbf{n} d \Gamma = - \int_{\Gamma} N^* (D \frac{\beta}{\alpha} N) \mathbf{n} d \Gamma U + \int_{\Gamma} N^* D \frac{c}{\alpha} \mathbf{n} d \Gamma
 * \f]
 * 
 */    
class RobinBC : IFemBC
{
public:
//    void set(TemplateFEMNodalFunction<double,double> *fem, int geo, int func);
};

}
