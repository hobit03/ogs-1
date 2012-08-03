/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CG.h
 *
 * Created on 2011-09-27 by Thomas Fischer
 */

#ifndef CG_H_
#define CG_H_

namespace MathLib {

// forward declaration
template <typename PF_TYPE, typename IDX_TYPE> class CRSMatrix;

unsigned CG(CRSMatrix<double,unsigned> const * mat, double const * const b,
        double* const x, double& eps, unsigned& nsteps);

#ifdef _OPENMP
unsigned CGParallel(CRSMatrix<double,unsigned> const * mat, double const * const b,
        double* const x, double& eps, unsigned& nsteps, unsigned num_threads = 1);
#endif

} // end namespace MathLib

#endif /* SOLVER_H_ */
