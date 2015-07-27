/*
 * EMOptimizer.h
 *
 *  Created on: Sep 9, 2013
 *      Author: afriesen
 */

#ifndef RDIS_EMOPTIMIZER_H_
#define RDIS_EMOPTIMIZER_H_

#include "common.h"
#include "RDISOptimizer.h"
#include "gmm/GMMFunction.h"


namespace rdis {

class EMOptimizer {
public:
	EMOptimizer( GMMFunction & fgmm_, Numeric epsilon_ = 1e-8 );
	~EMOptimizer();

	Numeric optimize( long long int maxiters = 10000, const State * xinitial = NULL,
			bool printdbg = true,
			Clock::time_point timeout = Clock::time_point::max() );

public:
	Numeric getOptVal() const { return fopt; }
	const State & getOptState() const { return xopt; }
	const Duration & getRuntime() const { return runtime; }
	size_t getNumIterations() const { return numIterations; }

private:
	void EStep();
	bool MStep();

private:
	// compute the log likelihood of the data given the current var assignments
	Numeric computeLL();

	// compute a_k * f(x_i | \mu_k, \sigma_k) -- assuming all vars are already
	// assigned and can be eval'd
	Numeric feval( size_t i, size_t k );

	// called from the M-step -- compute the maximum of the respective variables
	// given the results of the E-step
	bool maximizeWeights();
	bool maximizeMu();
	bool maximizeCov();

private:
	GMMFunction & fgmm;
	Numeric epsilon;

	// the set of posterior probabilities (indexed by varid i=1,...,n and
	// cluster id k=1,...,k)
	std::vector< std::vector< Numeric > > posteriorProbs;

	// optimal value computed
	Numeric fopt;

	// state at which the optimal value was computed
	State xopt;

	// total EM runtime
	Duration runtime;

	// total number of iterations taken
	size_t numIterations;
};

} // namespace rdis

#endif // RDIS_EMOPTIMIZER_H_
