/*
 * solver.h
 *
 *  Created on: Sep 27, 2011
 *      Author: TF
 */

#ifndef CG_H_
#define CG_H_

namespace MathLib {

// forward declaration
template <typename PF_TYPE, typename IDX_TYPE> class CRSMatrix;

unsigned CG(CRSMatrix<double,unsigned> const * mat, double const * const b,
		double* const x, double& eps, unsigned& nsteps, unsigned num_threads = 1);

} // end namespace MathLib

#endif /* SOLVER_H_ */
