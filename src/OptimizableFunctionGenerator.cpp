/*
 * OptimizableFunctionGenerator.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: afriesen
 */

#include "OptimizableFunctionGenerator.h"

#include "util/utility.h"
#include "VariableDomain.h"
#include "NonlinearProductFactor.h"
#include "PolynomialFunction.h"
#include "SimpleSumFactor.h"

#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/format.hpp>

#include BOOSTPATH/random/uniform_01.hpp>
#include BOOSTPATH/random/uniform_int_distribution.hpp>
#include BOOSTPATH/random/uniform_real_distribution.hpp>
#include BOOSTPATH/random/bernoulli_distribution.hpp>

#include BOOSTPATH/numeric/ublas/matrix.hpp>
#include BOOSTPATH/numeric/ublas/io.hpp>


namespace rdis {

BOOSTNS::random::mt19937 * OptimizableFunctionGenerator::randNumGen =
		new BOOSTNS::random::mt19937( 1351841847 );


OptimizableFunctionGenerator::OptimizableFunctionGenerator() {}
OptimizableFunctionGenerator::~OptimizableFunctionGenerator() {}


OptimizableFunctionSP OptimizableFunctionGenerator::generate(
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
		) {

	using namespace BOOSTNS;
	using namespace BOOSTNS::random;
	BOOSTNS::random::mt19937 & rng = *randNumGen;
	uniform_01<> unif01;

	assert( numVarsPerFact <= numVars );

	VariableDomain defaultDom( domain );
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( defaultDom );
	VariablePtrVec & vars = poly->variables;
	FactorPtrVec & factors = poly->factors;

	// create the empty factors
	for ( FactorID fid = 0; fid < (FactorID) numFactors; ++fid ) {
		Numeric fconst =
				( unif01( rng ) > factConstPct ) ? 0 : genConst( rng );
		Numeric fexp =
				( unif01( rng ) > factExpPct ) ? 1 : genExponent( rng );
		Numeric fcoeff =
				( unif01( rng ) > factCoeffPct ) ? 1 : genCoeff( rng );

		Factor * f = new SimpleSumFactor( fid, fconst, fexp, fcoeff );
		factors.push_back( f );
	}

	// variable naming scheme
	format fvn( "x%1%" );
	VariableID unique_vid( 0 );

	// place the backbone variables
	for ( VariableID vid = 0; vid < (long long int) numInBackbone; ++vid ) {
		string vname = str( fvn % vid );
		Variable * v = NULL;
		for ( FactorID fid = 0; fid < numFactors; ++fid ) {
			if ( factors[fid]->numVars() < numVarsPerFact &&
					unif01( rng ) < backbonePct ) {
				if ( v == NULL ) v = poly->addVariable( vname, "", unique_vid );
				Numeric exp = ( unif01( rng ) > varExponentPct ) ? 1 : genExponent( rng );
				Numeric coeff = ( unif01( rng ) > varCoeffPct ) ? 1 : genCoeff( rng );
				putVarInFactor( v, factors[fid], exp, coeff );
			}
		}
	}

	// create a list of the non-BB variables to avoid duplicates
	BOOSTNS::container::vector< VariableID > varIDs( numVars - numInBackbone );
	for ( VariableID i = numInBackbone; i < (VariableCount) numVars; ++i ) {
		varIDs[i - numInBackbone] = i;
		assert( i == unique_vid );
		poly->addVariable( str( fvn % i ), "", unique_vid );
	}

	// fill the remaining slots in each factor with randomly selected variables
	for ( Factor * f : factors ) {
		VariableCount needs = std::min( numVarsPerFact - f->numVars(), varIDs.size() );
		permute( varIDs, rng );
		for ( VariableCount i = 0; i < needs; ++i ) {
			// randomly choose unique variables from the non-backbone set
			VariableID vid = varIDs[i];
			Variable * v = vars[vid];
			assert( v != NULL && v->getID() == vid );

			// add it to this factor with randomly determined coeff & exponent
			Numeric exp = ( unif01( rng ) > varExponentPct ) ? 1 : genExponent( rng );
			Numeric coeff = ( unif01( rng ) > varCoeffPct ) ? 1 : genCoeff( rng );
			putVarInFactor( v, f, exp, coeff );
		}
	}

	// rename the used variables to align their names with their indices
	for ( Variable * v : poly->variables ) {
		v->setName( str( fvn % v->getID() ) );
	}

	poly->init();

	return poly;
}


void OptimizableFunctionGenerator::putVarInFactor( Variable * v,
		Factor * f, Numeric exp, Numeric coeff ) {

	// link this variable and factor
	dynamic_cast< SimpleSumFactor * >( f )->addVariable( v, exp, coeff );
}


Numeric OptimizableFunctionGenerator::genConst( BOOSTNS::random::mt19937 & rng ) {
	BOOSTNS::random::uniform_real_distribution<> unif( -2, 2 );
	return unif( rng );
}

Numeric OptimizableFunctionGenerator::genExponent( BOOSTNS::random::mt19937 & rng ) {
	BOOSTNS::random::uniform_int_distribution<> unif( 0, 5 );
	Numeric e( unif( rng ) );
	return e;
}

Numeric OptimizableFunctionGenerator::genCoeff( BOOSTNS::random::mt19937 & rng ) {
	BOOSTNS::random::uniform_real_distribution<> unif( -3, 3 );
	Numeric c( unif( rng ) );
	return c;
}


void OptimizableFunctionGenerator::createFactors( OptimizableFunctionSP & func,
		int nfactors ) {
	FactorPtrVec & factors = func->factors;
	for ( int i = 0; i < nfactors; ++i ) {
		factors.push_back( new SimpleSumFactor( i ) );
	}
}


void OptimizableFunctionGenerator::createVars( OptimizableFunctionSP & func,
		int nvars ) {
	BOOSTNS::format fvn( "x%1%" );
	VariableID unique_vid( 0 );
	for ( int i = 0; i < nvars; ++i ) {
		func->addVariable( BOOSTNS::str( fvn % unique_vid ), "", unique_vid );
	}
}


//#ifdef DEBUG
OptimizableFunctionSP OptimizableFunctionGenerator::makeSpecific() {

//	f(x,y) = -3.53738.[3.45169 + 2.45731.x^5 + 2.89465.y^3]^2
//	DIFF: max of 384719 at (2.74307, -2.61747)
//	RANGE: max of 392392 at ((2.74307, 1.74408), -2.61747) in
//		range:  [-1.05307, 2.74307], [-3.1599, -2.61747]

//	(-3)[ 5 + x0^(0) + (-1)x1^(3) ]^(0)
//		range [-0.5, 0.9]

	VariableDomain domain( "-3:4" );
	OptimizableFunctionSP poly(
			BOOSTNS::make_shared< PolynomialFunction >( domain ) );
	FactorPtrVec & factors = poly->factors;

	Factor * f = new SimpleSumFactor( 0, 7.92938, 2, 1.37471 );
	factors.push_back( f );

	// variable naming scheme
	BOOSTNS::format fvn( "x%1%" );
	VariableID unique_vid( 0 );
	string vname = BOOSTNS::str( fvn % unique_vid );
	Variable * v = poly->addVariable( vname, "", unique_vid );
	putVarInFactor( v, f, 5, 1.27546 );

	vname = BOOSTNS::str( fvn % unique_vid );
	v = poly->addVariable( vname, "", unique_vid );
	putVarInFactor( v, f, 1, 0.999404 );

	// initialize all of the factors and variables
	poly->init();

	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeSimplePoly() {

	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-3:4" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 2 );
	createVars( poly, 3 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0] );

	putVarInFactor( variables[0], factors[1] );
	putVarInFactor( variables[2], factors[1] );

	// initialize all of the factors and variables
	poly->init();

	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeMinStateWrongPoly() {

	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "1:4" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 1 );
	createVars( poly, 2 );

	variables[0]->setDomain( VariableDomain( -2, 1 ));
	variables[1]->setDomain( VariableDomain( 0.25, 0.5 ));

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0] );

	((SimpleSumFactor &) *factors[0]).setExponent( 2 );
	((SimpleSumFactor &) *factors[0]).setCoeff( -1 );

	// initialize all of the factors and variables
	poly->init();

	return poly;
}



OptimizableFunctionSP OptimizableFunctionGenerator::makeNonDecompPoly() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-2:4" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 3 );
	createVars( poly, 3 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0] );

	putVarInFactor( variables[0], factors[1] );
	putVarInFactor( variables[2], factors[1] );

	putVarInFactor( variables[1], factors[2] );
	putVarInFactor( variables[2], factors[2] );

	poly->init();

	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeSetToConstFactor() {
	OptimizableFunctionSP poly = makeNonDecompPoly();
	FactorPtrVec & factors = poly->factors;
	dynamic_cast< SimpleSumFactor * >( factors[2] )->setExponent( 0 );
	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeSetToConstFactor2() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "1:4" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 3 );
	createVars( poly, 3 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0] );

	putVarInFactor( variables[0], factors[1] );
	putVarInFactor( variables[2], factors[1] );

	putVarInFactor( variables[1], factors[2], -6 );
	putVarInFactor( variables[2], factors[2], -6 );

	// initialize all of the factors and variables
	poly->init();

	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeTreePoly() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "1:4" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 4 );
	createVars( poly, 5 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0] );

	putVarInFactor( variables[0], factors[1] );
	putVarInFactor( variables[2], factors[1] );

	putVarInFactor( variables[1], factors[2] );
	putVarInFactor( variables[3], factors[2] );

	putVarInFactor( variables[2], factors[3] );
	putVarInFactor( variables[4], factors[3] );

	poly->init();

	return poly;
}


//OptimizableFunctionSP OptimizableFunctionGenerator::makeTreePoly2() {
//	OptimizableFunctionSP poly =
//			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "1:4" ) );
//	FactorPtrVec & factors = poly->factors;
//	VariablePtrVec & variables = poly->variables;
//
//	createFactors( poly, 4 );
//	createVars( poly, 5 );
//
//	putVarInFactor( variables[0], factors[0] );
//	putVarInFactor( variables[1], factors[0] );
//
//	putVarInFactor( variables[1], factors[1], -6 );
//	putVarInFactor( variables[2], factors[1], -6 );
////	putVarInFactor( variables[1], factors[1], 1, -0.5 );
////	putVarInFactor( variables[2], factors[1], 1, -0.5 );
////	((SimpleSumFactor &) *factors[1]).setConstant( 5 );
////	((SimpleSumFactor &) *factors[1]).setExponent( -6 );
//
//
//	putVarInFactor( variables[2], factors[2], 1, -0.5 );
//	putVarInFactor( variables[3], factors[2], 1, -0.5 );
//	((SimpleSumFactor &) *factors[2]).setConstant( 5 );
//	((SimpleSumFactor &) *factors[2]).setExponent( -6 );
//
//	putVarInFactor( variables[3], factors[3] );
//	putVarInFactor( variables[4], factors[3] );
//
//	// initialize all of the factors and variables
//	poly->init();
//
//	return poly;
//}


OptimizableFunctionSP OptimizableFunctionGenerator::make2ComponentPoly() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "1:4" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 4 );
	createVars( poly, 6 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[2], factors[0] );

	putVarInFactor( variables[0], factors[1] );
	putVarInFactor( variables[3], factors[1] );

	putVarInFactor( variables[1], factors[2] );
	putVarInFactor( variables[4], factors[2] );

	putVarInFactor( variables[1], factors[3] );
	putVarInFactor( variables[5], factors[3] );

	// initialize all of the factors and variables
	poly->init();

	return poly;
}



OptimizableFunctionSP OptimizableFunctionGenerator::makeTreePoly3() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-2:2" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 14 );
	createVars( poly, 11 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0] );

	putVarInFactor( variables[0], factors[1] );
	putVarInFactor( variables[2], factors[1] );

	putVarInFactor( variables[0], factors[2] );
	putVarInFactor( variables[3], factors[2] );

	putVarInFactor( variables[0], factors[3] );
	putVarInFactor( variables[4], factors[3] );

	putVarInFactor( variables[1], factors[4] );
	putVarInFactor( variables[2], factors[4] );

	putVarInFactor( variables[1], factors[5] );
	putVarInFactor( variables[3], factors[5] );

	putVarInFactor( variables[1], factors[6] );
	putVarInFactor( variables[4], factors[6] );

	putVarInFactor( variables[2], factors[7] );
	putVarInFactor( variables[5], factors[7] );

	putVarInFactor( variables[2], factors[8] );
	putVarInFactor( variables[6], factors[8] );

	putVarInFactor( variables[3], factors[9] );
	putVarInFactor( variables[7], factors[9] );

	putVarInFactor( variables[3], factors[10] );
	putVarInFactor( variables[8], factors[10] );

	putVarInFactor( variables[4], factors[11] );
	putVarInFactor( variables[9], factors[11] );

	putVarInFactor( variables[4], factors[12] );
	putVarInFactor( variables[10], factors[12] );

	putVarInFactor( variables[0], factors[13] );
	putVarInFactor( variables[5], factors[13] );

	for ( int i = 0; i < 14; ++i ) {
		dynamic_cast< SimpleSumFactor * >( factors[i] )->setExponent( 2 );
	}

	poly->init();
	return poly;
}



OptimizableFunctionSP OptimizableFunctionGenerator::makeCrossPoly() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "1:1.3" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 4 );
	createVars( poly, 8 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[2], factors[0] );
	putVarInFactor( variables[4], factors[0] );
	putVarInFactor( variables[6], factors[0] );

	putVarInFactor( variables[1], factors[1] );
	putVarInFactor( variables[3], factors[1] );
	putVarInFactor( variables[5], factors[1] );
	putVarInFactor( variables[7], factors[1] );

	putVarInFactor( variables[0], factors[2] );
	putVarInFactor( variables[2], factors[2] );
	putVarInFactor( variables[5], factors[2] );
	putVarInFactor( variables[7], factors[2] );

	putVarInFactor( variables[1], factors[3] );
	putVarInFactor( variables[3], factors[3] );
	putVarInFactor( variables[4], factors[3] );
	putVarInFactor( variables[6], factors[3] );


	// initialize all of the factors and variables
	poly->init();
	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makePowellsFunction() {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-5:5" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 4 );
	createVars( poly, 4 );

	putVarInFactor( variables[0], factors[0] );
	putVarInFactor( variables[1], factors[0], 1, 10 );
	dynamic_cast< SimpleSumFactor * >( factors[0] )->setExponent( 2 );
	dynamic_cast< SimpleSumFactor * >( factors[0] )->setCoeff( 0.5 );

	putVarInFactor( variables[2], factors[1] );
	putVarInFactor( variables[3], factors[1], 1, -1 );
	dynamic_cast< SimpleSumFactor * >( factors[1] )->setExponent( 2 );
	dynamic_cast< SimpleSumFactor * >( factors[1] )->setCoeff( 5.0 * 0.5 );

	putVarInFactor( variables[1], factors[2] );
	putVarInFactor( variables[2], factors[2], 1, -2 );
	dynamic_cast< SimpleSumFactor * >( factors[2] )->setExponent( 4 );
	dynamic_cast< SimpleSumFactor * >( factors[2] )->setCoeff( 0.5 );

	putVarInFactor( variables[0], factors[3] );
	putVarInFactor( variables[3], factors[3], 1, -1 );
	dynamic_cast< SimpleSumFactor * >( factors[3] )->setExponent( 4 );
	dynamic_cast< SimpleSumFactor * >( factors[3] )->setCoeff( 10.0 * 0.5 );


	// initialize all of the factors and variables
	poly->init();
	return poly;
}


//OptimizableFunctionSP OptimizableFunctionGenerator::makeColville() {
////	y  = 100*(x(1)^2-x(2))^2 + (x(1)-1)^2 + (x(3)-1)^2 + 90*(x(3)^2-x(4))^2 + ...
////		10.1*((x(2)-1)^2+(x(4)-1)^2)+19.8*(x(2)^-1)*(x(4)-1);
//
////	y  =
////	100 * (x1^2-x2)^2
////	+ (x1-1)^2
////	+ (x3-1)^2
////	+ 90*(x3^2-x4)^2
////	+ 10.1*((x2-1)^2 + (x4-1)^2)  --> 10.1*(x2-1)^2 + 10.1*(x4-1)^2
////	+ 19.8*(x2^-1)*(x4-1) -->
//	OptimizableFunctionSP poly =
//			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-10:10" ) );
//	FactorPtrVec & factors = poly->factors;
//	VariablePtrVec & variables = poly->variables;
//
//	createFactors( poly, 7 );
//	createVars( poly, 4 );
//
//	//	100 * (x1^2-x2)^2
//	putVarInFactor( variables[0], factors[0], 2 );
//	putVarInFactor( variables[1], factors[0], 1, -1 );
//	dynamic_cast< SimpleSumFactor * >( factors[0] )->setExponent( 2 );
//	dynamic_cast< SimpleSumFactor * >( factors[0] )->setCoeff( 100 );
//
//	//	+ (x1-1)^2
//	putVarInFactor( variables[0], factors[1] );
//	dynamic_cast< SimpleSumFactor * >( factors[1] )->setExponent( 2 );
//	dynamic_cast< SimpleSumFactor * >( factors[1] )->setConstant( -1 );
//
//	//	+ (x3-1)^2
//	putVarInFactor( variables[2], factors[2] );
//	dynamic_cast< SimpleSumFactor * >( factors[2] )->setExponent( 2 );
//	dynamic_cast< SimpleSumFactor * >( factors[2] )->setConstant( -1 );
//
//	//	+ 90*(x3^2-x4)^2
//	putVarInFactor( variables[2], factors[3], 2 );
//	putVarInFactor( variables[3], factors[3], 1, -1 );
//	dynamic_cast< SimpleSumFactor * >( factors[3] )->setExponent( 2 );
//	dynamic_cast< SimpleSumFactor * >( factors[3] )->setCoeff( 90 );
//
//	//	+ 10.1*((x2-1)^2 + (x4-1)^2)
//	putVarInFactor( variables[1], factors[4] );
//	dynamic_cast< SimpleSumFactor * >( factors[4] )->setExponent( 2 );
//	dynamic_cast< SimpleSumFactor * >( factors[4] )->setConstant( -1 );
//	dynamic_cast< SimpleSumFactor * >( factors[4] )->setCoeff( 10.1 );
//
//	putVarInFactor( variables[3], factors[5] );
//	dynamic_cast< SimpleSumFactor * >( factors[5] )->setExponent( 2 );
//	dynamic_cast< SimpleSumFactor * >( factors[5] )->setConstant( -1 );
//	dynamic_cast< SimpleSumFactor * >( factors[5] )->setCoeff( 10.1 );
//
//	//	+ 19.8*(x2^-1)*(x4-1)
//	NonlinearProductFactor * spf =
//			new NonlinearProductFactor( factors[6]->getID(), 19.8, false );
//	delete factors[6];
//	factors[6] = spf;
//	spf->addVariable( variables[1], -1 );
//	spf->addVariable( variables[3], 1, -1 );
//
//	poly->init();
//
//	return poly;
//}


OptimizableFunctionSP OptimizableFunctionGenerator::makeRosenbrock( const int N ) {
//	J=zeros(m,n);
//
//	for i=1:m/2
//
//	   if (option==1 | option==3)
//	        fvec(2*i-1)=10*(x(2*i)-x(2*i-1)^2);
//	        fvec(2*i)=1-x(2*i-1);
//	    else fvec='?';
//	   end;
//
//	   if (option==2 | option==3)
//	        J(2*i-1,2*i-1) = -20*x(2*i-1);
//	        J(2*i-1,2*i)   = 10;
//	        J(2*i,2*i-1)   = -1;
//	    else J='?';
//	   end;
//
//	end;

	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-5:5" ) );
	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	createFactors( poly, 2*(N-1) );
	createVars( poly, N );

	int fid = 0;
	for ( int i = 0; i < N-1; ++i ) {
		Factor * f1 = factors[fid++];
		Factor * f2 = factors[fid++];

		putVarInFactor( variables[i], f1, 1, -1 );
		dynamic_cast< SimpleSumFactor * >( f1 )->setExponent( 2 );
		dynamic_cast< SimpleSumFactor * >( f1 )->setConstant( 1 );

		putVarInFactor( variables[i+1], f2 );
		putVarInFactor( variables[i], f2, 2, -1 );
		dynamic_cast< SimpleSumFactor * >( f2 )->setExponent( 2 );
		dynamic_cast< SimpleSumFactor * >( f2 )->setCoeff( 100 );
	}

	poly->init();
	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeHighDimSinusoid(
		const VariableCount treeHeight, const VariableCount branches,
		VariableCount maxArity, const bool allowOddArityFactors ) {

	const Numeric twopi = 2.000001*3.141592653;
	const NumericInterval samplingInterval( -twopi, twopi );

	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >(
					VariableDomain( BOOSTNS::str(
							BOOSTNS::format( "-%1%:%1%" ) % ( 10*twopi ) ) ) );

	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	maxArity = std::min( maxArity, treeHeight+1 );

	const VariableCount h = treeHeight;
	const VariableCount k = branches;

	// formula for complete k-ary tree
	const VariableCount nvars = ( k == 1 ? h+1 :
			( round( std::pow( (Numeric) k, (Numeric) h+1 ) ) - 1 ) / ( k - 1 ) );

	std::cout << "tree: h = " << h << ", k = " << k << ", nv = " << nvars << std::endl;

	createVars( poly, nvars );

	// set the sampling interval for this var
	for ( Variable * v : variables ) {
		v->setSamplingInterval( samplingInterval );
	}

	FactorID fid = 0;
	VariablePtrVec tmpvars;

	for ( VariableCount ar = 1; ar <= maxArity; ++ar ) {
		if ( ar > 1 && ( ar & 0x01 ) > 0 && !allowOddArityFactors ) continue;

//		std::cout << "arity " << ar << std::endl;
		VariableCount lasth = h;

		for ( VariableID vid = nvars-1; vid >= 0; --vid ) {
			VariableID lastVidAtNextH = k == 1 ? lasth-1 :
					( round( std::pow( (Numeric) k, (Numeric) lasth ) )-1.0 ) /
						( k-1.0 ) - 1.0;
			VariableCount varheight = ( vid > lastVidAtNextH ? lasth : --lasth );

//			std::cout << "var " << vid << " at height " << varheight <<
//					", factor " << fid << " (arity " << ar << ") -- " <<
//					lastVidAtNextH << std::endl;

			if ( varheight+1 < ar ) continue;

			tmpvars.clear();
			VariableID cur = vid, maxvid = vid;

			// walk up the tree and get the ancestors up to the current arity
			for ( VariableCount c = 0; c < ar; ++c ) {
				assert( cur >= 0 );
				tmpvars.push_back( variables[cur] );
				if ( cur > maxvid ) maxvid = cur;
				VariableID par = std::floor( ( (Numeric) cur - 1.0 ) / (Numeric) k );
//				std::cout << "added var " << cur << " (next " << par << ")" << std::endl;
				cur = par;
			}

			Numeric coeff = ar > 1 ? 12 : 0.6;
			NonlinearProductFactor *nlpf =
					new NonlinearProductFactor( fid++, coeff, false, maxvid+1 );
			factors.push_back( nlpf );

			// add the variables
			while ( !tmpvars.empty() ) {
				nlpf->addVariable( tmpvars.back(), 1, 0, ar > 1 );
				tmpvars.pop_back();
			}

			assert( (VariableCount) nlpf->getVariables().size() == ar );
		}
	}

	for ( VariableID vid = 0; vid < nvars; ++vid ) {
		NonlinearProductFactor *nlpf = new NonlinearProductFactor( fid++, 0.1, false );
		factors.push_back( nlpf );
		nlpf->addVariable( variables[vid], 2, 0, false );
	}

	std::cout << "created " << factors.size() << " factors" << std::endl;

//	for ( Factor * f : factors ) {
//		std::cout << "factor " << f->getID() << " has " << f->getVariables().size()
//				<< " vars: ";
//		for ( Variable * v : f->getVariables() ) std::cout << v->getID() << ", ";
//		std::cout << std::endl;
//	}

	poly->init();

	return poly;
}


OptimizableFunctionSP OptimizableFunctionGenerator::makeRandomConvex(
		const VariableCount numVars, const Numeric percentSparse, State & xinit ) {
	OptimizableFunctionSP poly =
			BOOSTNS::make_shared< PolynomialFunction >( VariableDomain( "-2:2" ) );

	createVars( poly, numVars );

	FactorPtrVec & factors = poly->factors;
	VariablePtrVec & variables = poly->variables;

	FactorID fid( 0 );
	BOOSTNS::numeric::ublas::matrix< Numeric > M( numVars, numVars );
	getSparsePSDMatrix( M, percentSparse );

	BOOSTNS::numeric::ublas::zero_vector< Numeric > b( numVars );

	for ( VariableCount i = 0; i < (VariableCount) M.size1(); ++i ) {
		for ( VariableCount j = 0; j < (VariableCount) M.size2(); ++j ) {

			if ( i > j || M( i, j ) == 0 ) continue;

			Numeric coeff( M( i, j ) * ( i == j ? 1 : 2 ) );
			factors.push_back( new NonlinearProductFactor( fid++, coeff ) );
			NonlinearProductFactor & nlpf( dynamic_cast< NonlinearProductFactor & >(
						*factors.back() ) );
			if ( i == j ) {
				nlpf.addVariable( variables[i], 2, b( i ) );
			} else {
				nlpf.addVariable( variables[i], 1, b( i ) );
				nlpf.addVariable( variables[j], 1, b( j ) );
			}
		}
	}

	using namespace BOOSTNS::random;
	xinit.resize( numVars, 0 );
	for ( VariableCount i = 0; i < numVars; ++i ) {
		uniform_real_distribution<> unif( variables[i]->getDomain().min(),
										  variables[i]->getDomain().max() );
		xinit[i] = unif( *randNumGen );
	}

	return poly;
}


void OptimizableFunctionGenerator::getSparsePSDMatrix(
		BOOSTNS::numeric::ublas::matrix< Numeric > & M, Numeric percentSparse ) {

	using namespace BOOSTNS::random;
	using namespace BOOSTNS::numeric::ublas;

	const Numeric scale( 1.0 );
	uniform_real_distribution<> unif( -scale, scale );
	bernoulli_distribution<> flip( percentSparse );
//	bernoulli_distribution<> flip( 1.0 - percentSparse );
	mt19937 & rng = *randNumGen;

	const VariableCount rows( M.size1() ), cols( M.size2() );

	for ( VariableCount i = 0; i < rows; ++i ) {
		for ( VariableCount j = 0; j < cols; ++j ) {

			if ( i < j ) {
				M( i, j ) = 0;

			} else if ( i == j ) {
				M( i, j ) = unif( rng ) + 2*scale;

			} else {
				// otherwise, off-diagonal element in lower triangular part
				M( i, j ) = unif( rng );

//				// flip a coin to determine how to set this value
//				if ( flip( rng ) ) M( i, j ) = unif( rng );
//				else M( i, j ) = 0;
			}
		}
	}

//	matrix< Numeric > Mt( trans( M ) );
	M = prod( M, trans( M ) );

	// "sparsify" the matrix
	for ( VariableCount i = 0; i < rows; ++i ) {
		for ( VariableCount j = 0; j < cols; ++j ) {
			if ( i < j || i == j ) continue;
			if ( flip( rng ) ) M( i, j ) = M( j, i ) = 0;
		}
	}
}


//#endif // DEBUG

} // namespace rdis
