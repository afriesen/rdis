/*
 * LMSubspaceOptimizer.h
 *
 *  Created on: Apr 30, 2014
 *      Author: afriesen
 */

#ifndef RDIS_LMSUBSPACEOPTIMIZER_H_
#define RDIS_LMSUBSPACEOPTIMIZER_H_

#include "common.h"
#include "SubspaceOptimizer.h"

#include "OptimizableFunction.h"

namespace rdis {

// Levenberg-Marquardt based subspace optimizer
class LMSubspaceOptimizer
		: public rdis::SubspaceOptimizer {

	typedef rdis::SubspaceOptimizer super;

public:
	LMSubspaceOptimizer( OptimizableFunction & f_ );
	virtual ~LMSubspaceOptimizer();

public:
	// parameters are the vars and factors that define this subfunction, the
	// starting values of those variables in the same order as the <vars>
	// vector (will contain the final values after optimization completes), and
	// the change in f(x) as a result of the optimization = f(x_end) - f(x_init)
	virtual Numeric optimize( const VariablePtrVec & vars,
			const FactorPtrVec & factors, NumericVec & xval,
			Numeric & deltaFval, const bool printdbg );

protected:

	// store for memory efficiency
//	NumericVec initstate;
	NumericVec initeval;
	NumericVec finaleval;
	PartialGradient pgtemp;
};

namespace LMSSOpt {

// struct containing auxiliary data for function evaluation
struct AuxData {
	AuxData( const VariablePtrVec & va, const FactorPtrVec & fac,
			OptimizableFunction & func,
			PartialGradient & pgt_ )
		: vars( va )
		, factors( fac )
		, f( func )
		, pgt( pgt_ )
	{}

	const VariablePtrVec & vars;
	const FactorPtrVec & factors;
	OptimizableFunction & f;

	// store for memory efficiency
	PartialGradient & pgt;
};


// evaluate the function and the jacobian of the function at the specified
// point p
void evalFunc( double * p, double * hx, int m, int n, void * adata );
void evalJacf( double * p, double * j, int m, int n, void * adata );

void quickAssignVals( AuxData & aux, int m, double * p );

} // namespace LMSSOpt

} // namespace rdis
#endif // RDIS_LMSUBSPACEOPTIMIZER_H_
