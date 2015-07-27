/*
 * EMSubspaceOptimizer.h
 *
 *  Created on: May 8, 2014
 *      Author: afriesen
 */

#ifndef EMSUBSPACEOPTIMIZER_H_
#define EMSUBSPACEOPTIMIZER_H_

#include "SubspaceOptimizer.h"
#include "OptimizableFunction.h"

namespace rdis {

// runs expectation-maximization (EM) algorithm on the specified subspace to
// optimize it
class EMSubspaceOptimizer
		: public SubspaceOptimizer {

	typedef SubspaceOptimizer super;

public:
	EMSubspaceOptimizer( OptimizableFunction & f_, size_t maxiters_ = 50,
			Numeric ftol_ = 1e-10, bool sanitizeVals = true,
			bool useFactorCaching_ = true );
	virtual ~EMSubspaceOptimizer();

public:
	virtual Numeric optimize( const VariablePtrVec & vars,
			const FactorPtrVec & factors, NumericVec & xinit,
			Numeric & deltaFval, const bool printdbg );

protected:
	// compute the log likelihood of this model
	Numeric computeLL( const FactorPtrVec & factors );

protected:
	// make sure that variable values are always within the variable's domain
	bool sanitizeVals;

	// true to not use factor caching (tradeoff between not notifying vars when
	// assigning them and having to re-evaluate factors each time)
	bool useFactorCaching;

//	// temporary vector and gradient storage to reduce memory overhead
//	IntervalVec ivtmp;
//	PartialGradient pgtmp;
//	PartialGradient pgtmp2;

};

} // namespace rdis
#endif // EMSUBSPACEOPTIMIZER_H_
