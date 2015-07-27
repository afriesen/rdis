/*
 * OptimizableFunctionGenerator.h
 *
 *  Created on: Jan 29, 2013
 *      Author: afriesen
 */

#ifndef RDIS_OPTIMIZABLEFUNCTIONGENERATOR_H_
#define RDIS_OPTIMIZABLEFUNCTIONGENERATOR_H_

#include "common.h"
#include "OptimizableFunction.h"

#include BOOSTPATH/random/mersenne_twister.hpp>
#include BOOSTPATH/numeric/ublas/fwd.hpp>

namespace rdis {

class OptimizableFunctionGenerator {
public:
	OptimizableFunctionGenerator();
	~OptimizableFunctionGenerator();

	static OptimizableFunctionSP generate(
			string domain, 			// domain of the variables
			size_t numFactors, 		// number of factors in this polynomial
			size_t numVars,			// number of variables
			size_t numVarsPerFact,	// number of variables in each factor
			size_t numInBackbone, 	// number of variables in the "backbone"
			float backbonePct,		// the percent of factors that BB vars are in
			float varExponentPct,	// the percent of vars with exponents
			float varCoeffPct,		// the percent of vars with coefficients
			float factConstPct,		// the percent of factors with constants
			float factExpPct,		// the percent of factors with exponents
			float factCoeffPct		// the percent of factors with coefficients
		);

	static OptimizableFunctionSP makeSpecific();
	static OptimizableFunctionSP makeSimplePoly();
	static OptimizableFunctionSP makeMinStateWrongPoly();
	static OptimizableFunctionSP makeNonDecompPoly();
	static OptimizableFunctionSP makeSetToConstFactor();
	static OptimizableFunctionSP makeSetToConstFactor2();
	static OptimizableFunctionSP make2ComponentPoly();
	static OptimizableFunctionSP makeTreePoly();
//	static OptimizableFunctionSP makeTreePoly2();
	static OptimizableFunctionSP makeTreePoly3();
	static OptimizableFunctionSP makeCrossPoly();
	static OptimizableFunctionSP makePowellsFunction();
	static OptimizableFunctionSP makeRosenbrock( const int N );

	static OptimizableFunctionSP makeHighDimSinusoid( const VariableCount nvars,
													  const VariableCount height, VariableCount maxArity,
													  const bool oddArityFactors );

	static OptimizableFunctionSP makeRandomConvex( VariableCount numVars,
			Numeric percentSparse, State & xinit );


private:
	static void putVarInFactor( Variable * v, Factor * f,
			Numeric exp = 1, Numeric coeff = 1 );

	static Numeric genConst( BOOSTNS::random::mt19937 & rng );
	static Numeric genExponent( BOOSTNS::random::mt19937 & rng );
	static Numeric genCoeff( BOOSTNS::random::mt19937 & rng );

	static void createFactors( OptimizableFunctionSP & func, int nfactors );
	static void createVars( OptimizableFunctionSP & func, int nvars );

	static void getSparsePSDMatrix( BOOSTNS::numeric::ublas::matrix< Numeric > & M,
			Numeric percentSparse );

	static BOOSTNS::random::mt19937 * randNumGen;
};

} // namespace rdis

#endif // RDIS_OPTIMIZABLEFUNCTIONGENERATOR_H_
