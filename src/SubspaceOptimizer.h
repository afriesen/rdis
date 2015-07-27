/*
 * SubspaceOptimizer.h
 *
 *  Created on: Mar 25, 2014
 *      Author: afriesen
 */

#ifndef RDIS_SUBSPACEOPTIMIZER_H_
#define RDIS_SUBSPACEOPTIMIZER_H_

#include "common.h"

#include "Variable.h"
#include "Factor.h"
#include "OptimizableFunction.h"

#include BOOSTPATH/program_options.hpp>

namespace rdis {

// a simple base for classes that perform a basic optimization procedure over a
// subset of the variables and a subset of the factors/terms of the main
// objective function being optimized -- typically used to choose values
class SubspaceOptimizer {
public:
	SubspaceOptimizer( OptimizableFunction & f_ );
	virtual ~SubspaceOptimizer();

public:
	virtual void setParameters(
			const BOOSTNS::program_options::variables_map & /*options*/ );

	// parameters are the vars and factors that define this subfunction, the
	// starting values of those variables in the same order as the <vars>
	// vector (will contain the final values after optimization completes), and
	// the change in f(x) as a result of the optimization = f(x_end) - f(x_init)
	virtual Numeric optimize( const VariablePtrVec & vars,
			const FactorPtrVec & factors, NumericVec & xval,
			Numeric & deltaFval, const bool printdbg ) = 0;

protected:
	void quickAssignVals( const VariablePtrVec & vars,
			const NumericVec & xval, bool sanitizeVals );

protected:
	// this class is optimizing a sub-function of this function
	OptimizableFunction & f;

	// based on the semiring, negate evals and derivatives and do ascent
	const bool doAscent;

	// maximum number of iterations for each call to this optimizer
	size_t maxiters;

	// the tolerance in f() below which we assume convergence
	Numeric ftol;

};


} // namespace rdis

#endif // RDIS_SUBSPACEOPTIMIZER_H_
