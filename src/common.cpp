/*
 * common.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: afriesen
 */


#include "common.h"

namespace rdis {

const Numeric NumericMAX =  std::numeric_limits< Numeric >::max();
const Numeric NumericMIN = -std::numeric_limits< Numeric >::max();

//const Numeric NumericLARGE = 1e100;

namespace instrumentation {

size_t g_numFactEvals = 0;
size_t g_numVarOps = 0;

size_t g_numProducts = 0;
size_t g_numSums = 0;
//size_t g_numExp = 0;
//size_t g_numMax = 0;
Duration g_runTime = Duration::zero();

size_t g_numDownhillSteps = 0;
size_t g_numSteps = 0;

size_t g_cacheInserts = 0;
size_t g_cacheHits = 0;
size_t g_cacheLookups = 0;
size_t g_cacheMidInserts = 0;

size_t g_counter1 = 0;
size_t g_counter2 = 0;

void clearInstrumentation() {
	g_numFactEvals = 0;
	g_numVarOps = 0;
	g_numProducts = 0;
	g_numSums = 0;
//	g_numExp = 0;
//	g_numMax = 0;
	g_runTime = Duration::zero();

	g_numDownhillSteps = 0;
	g_numSteps = 0;

	g_cacheInserts = 0;
	g_cacheHits = 0;
	g_cacheLookups = 0;
	g_cacheMidInserts = 0;

	g_counter1 = 0;
	g_counter2 = 0;
}

} // namespace instrumentation
} // namespace rdis
