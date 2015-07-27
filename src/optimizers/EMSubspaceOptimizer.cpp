/*
 * EMSubspaceOptimizer.cpp
 *
 *  Created on: May 8, 2014
 *      Author: afriesen
 */

#include "optimizers/EMSubspaceOptimizer.h"

namespace rdis {

EMSubspaceOptimizer::EMSubspaceOptimizer( OptimizableFunction & f_,
		size_t maxiters_, Numeric ftol_, bool sanitizeVals_,
		bool useFactorCaching_ )
	: super( f_ )
	, sanitizeVals( sanitizeVals_ )
	, useFactorCaching( useFactorCaching_ )
{
	assert( !doAscent );
}

EMSubspaceOptimizer::~EMSubspaceOptimizer() {}


Numeric EMSubspaceOptimizer::optimize( const VariablePtrVec & vars,
		const FactorPtrVec & facs, NumericVec & xval, Numeric & deltaFval,
		const bool printdbg ) {

	bool result = true;
	Numeric curval( 0 ), prevval( std::numeric_limits< Numeric >::max() );

	const Clock::time_point starttime( Clock::now() );

	quickAssignVals( vars, xval, sanitizeVals );

	Numeric ferr = 0.0;
	Numeric initval = f.evalFactors( facs, ferr );
	prevval = curval = initval;

	VariableCount iter = 0;

	for ( ; ; ) {
//		if ( printdbg && iter == 0 ) {
//			std::clog << "iteration " << iter << ": log likelihood: " <<
//					curval << std::endl;
//		}

		// take an EM step (changes the var assignments)
		result = f.EMStep( vars, facs, sanitizeVals, useFactorCaching );

		// compute the new log likelihood
		prevval = curval;
		curval = f.evalFactors( facs, ferr, useFactorCaching );

		if ( !result ) {
			std::cout << "EM subspace opt. failed on iteration " << iter <<
					", log likelihood: " << curval << std::endl;
			break;
		}

		// if we've converged, get the final state and unassign all vars
		if ( ++iter >= (VariableCount) maxiters && maxiters > 0 ) break;
//		if ( !approxleq( curval, prevval, ftol ) ) break;
		if ( ( prevval - curval ) < ftol ) break;

//		break;
	}

	for ( size_t i = 0; i < vars.size(); ++i ) {
		xval[i] = vars[i]->eval();
	}

//	std::cout << "prevval " << prevval << ", curval " << curval << ", diff " << curval-prevval << std::endl;

	deltaFval = curval - initval;
	Duration runtime = ( Clock::now() - starttime );

	if ( printdbg ) {
	std::cout << "EM ss opt. (iters "<< iter << "): " << curval <<
			" (init: " << initval << ", diff: " << deltaFval << ") in " <<
			runtime.count() << " seconds" << std::endl;
	}

	return curval;
}

} // namespace rdis
