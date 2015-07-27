/*
 * BCDOptimizer.h
 *
 *  Created on: Jan 26, 2014
 *      Author: afriesen
 */

#ifndef RDIS_BCDOPTIMIZER_H_
#define RDIS_BCDOPTIMIZER_H_

#include "common.h"

#include "SubspaceOptimizer.h"

#include "external/include/minimize_nrc.h"

#include BOOSTPATH/random.hpp>

namespace rdis {

class BCDOptimizer {
public:
	BCDOptimizer( OptimizableFunction & f_, SubspaceOptimizer & ssopt_,
			size_t blocksize_, Numeric epsilon_ = 1e-8 );
	~BCDOptimizer();

	Numeric optimize( long long int maxitersCG = 100,
			long long int maxitersTotal = 10000, const State * xinitial = NULL,
			bool printdbg = true,
			Clock::time_point timeout = Clock::time_point::max(),
			size_t nRandomRestarts = 0,
			Clock::time_point actualStartTime = Clock::time_point::max() );

public:
	Numeric getOptVal() const { return fopt; }
	const State & getOptState() const { return xopt; }
	const Duration & getRuntime() const { return runtime; }
	size_t getNumIterations() const { return numIterations; }

	bool getTimedOut() const { return timedOut; }

private:
	void assignVars( const State & x );


private:
	OptimizableFunction & fnc;
	SubspaceOptimizer & ssopt;
	size_t blocksize;
	Numeric epsilon;

	bool timedOut;

	// optimal value computed
	Numeric fopt;

	// state at which the optimal value was computed
	State xopt;

	// total runtime
	Duration runtime;

	// total number of iterations taken
	size_t numIterations;

	BOOSTNS::random::mt19937 rng;

	std::vector< VariablePtrVec > blocks;
	std::vector< FactorPtrVec > fblocks;
};

} // namespace rdis

#endif // RDIS_BCDOPTIMIZER_H_
