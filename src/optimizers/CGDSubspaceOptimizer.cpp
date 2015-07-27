/*
 * CGDSubspaceOptimizer.cpp
 *
 *  Created on: Mar 25, 2014
 *      Author: afriesen
 */

#include "optimizers/CGDSubspaceOptimizer.h"

namespace rdis {

CGDSubspaceOptimizer::CGDSubspaceOptimizer( OptimizableFunction & f_ )
	: super( f_ )
{}

CGDSubspaceOptimizer::~CGDSubspaceOptimizer() {}


Numeric CGDSubspaceOptimizer::optimize( const VariablePtrVec & vars,
		const FactorPtrVec & gdfs, NumericVec & xval, Numeric & deltaFval,
		const bool printdbg ) {

	assert( xval.size() == vars.size() );

	// if nothing to optimize, return 0 and leave xval the same as it was
	if ( gdfs.empty() ) {
		deltaFval = 0;
		return 0;
	}

	// make sure the vars are already assigned the starting value
//	assert( vars[0]->isAssigned() && vars[0]->eval() == xval[0] );
	quickAssignVals( vars, xval, true );

	SubfunctionFD sfd( f, doAscent, vars, gdfs, ivtmp, pgtmp, pgtmp2 );

	const Numeric initialFval = sfd.coeff * sfd( xval );
	nvtmp.assign( xval.begin(), xval.end() );
	const NumericVec & initxval( nvtmp );
	nrc::Frprmn< SubfunctionFD > gdmin( sfd, maxiters, ftol );

	try {
		gdmin.minimize( xval );
	} catch ( const std::exception & e ) {
		if ( printdbg ) {
			std::cout << "caught exception thrown by NRC Frprmn: " <<
					e.what() << std::endl;
		}
	} catch ( const char * s ) {
		if ( printdbg ) {
			std::cout << "caught exception thrown by NRC Frprmn: " << s << std::endl;
		}

//		NumericVec deriv;
//		sfd.df( gdmin.p, deriv );
//		std::cout << "position: " << gdmin.p << std::endl;
//		std::cout << "derivative: " << deriv << std::endl;
	}

	// assign the vars their final vals (with sanitization)
	sfd.quickAssignVals( gdmin.p );

	Numeric fret = sfd.coeff * gdmin.fret;

	// if negative progress was made, restore the initial values
	if ( fret > initialFval ) {
		if ( printdbg ) {
			std::cout << "CGD SS: negative progress was made, so restoring initial values ("
					<< fret << " > " << initialFval << ")" << std::endl;
		}
//		fret = initialFval;
		sfd.quickAssignVals( initxval );
		fret = sfd.coeff * sfd( initxval );
		if ( !approxeq( fret, initialFval, 1e-10 ) ) {
			std::cout << "fval for initxvals does not equal initialFval: " <<
					fret << " != " << initialFval << " (diff: " <<
					fret - initialFval << ")" << std::endl;
		}
		assert( approxeq( fret, initialFval, 1e-10 ) );
	}

	// sanitize values to be within the domain of this variable (note that they
	// were sanitized when assigned, so just copy them out)
	for ( size_t i = 0; i < vars.size(); ++i ) {
		xval[i] = vars[i]->eval();
	}

//	nsteps = gdmin.iter;
	deltaFval = ( fret - initialFval );

	if ( printdbg ) {
		std::cout << "CGD subspace result (steps " << gdmin.iter << "): " <<
				fret << ", diff: " << deltaFval << ", init: " <<
				initialFval << std::endl;
	}

	return fret;
}



CGDSubspaceOptimizer::SubfunctionFD::SubfunctionFD(
		OptimizableFunction & func_, bool doAscent,
		const VariablePtrVec & vars_, const FactorPtrVec & facs_,
		IntervalVec & ivtemp, PartialGradient & pgtemp1_,
		PartialGradient & pgtemp2_ )
	: func( func_ )
	, coeff( doAscent ? -1 : 1 )
	, vars( vars_ )
	, facs( facs_ )
	, bounds( ivtemp )
	, pgtemp1( pgtemp1_ )
	, pgtemp2( pgtemp2_ ) {

	bounds.clear();
	for ( size_t i = 0; i < vars.size(); ++i ) {
		bounds.push_back( vars[i]->getDomain().interval() );
		// relax this! -- what is the input value? (assert var is assigned, and use that val...)
		assert( vars[i]->getDomain().subintervals().size() == 1 );
	}
}


Numeric CGDSubspaceOptimizer::SubfunctionFD::operator()( nrc::VecDoub_I & x ) {

	quickAssignVals( x );

	Numeric ferr = 0.0;
	Numeric result = func.evalFactors( facs, ferr, true );

	return coeff * result;
}


void CGDSubspaceOptimizer::SubfunctionFD::df( nrc::VecDoub_I & x,
		nrc::VecDoub_O & deriv ) {

	deriv.assign( vars.size(), 0 );

	quickAssignVals( x );

	// compute gradient
	pgtemp1.clear();
	func.computeGradient( facs, pgtemp1, false );

//	std::cout << "derivative of f() at " << x << " -- pg: ";
//	for ( auto it = pgtemp1.begin(); it != pgtemp1.end(); ++it ) {
//		std::cout << "(" << it->first << ", " << it->second << "), ";
//	}
//	std::cout << std::endl;

	// copy gradient into vector
	for ( size_t i = 0; i < vars.size(); ++i ) {
		auto it = pgtemp1.find( vars[i]->getID() );
		deriv[i] = coeff * ( it == pgtemp1.end() ? 0 : it->second );
	}
}


bool CGDSubspaceOptimizer::SubfunctionFD::quickAssignVals(
		const NumericVec & xval ) {

	bool res = true;

	for ( size_t i = 0; i < vars.size(); ++i ) {
		Numeric val = vars[i]->getDomain().closestVal( xval[i] );
		// TODO: Fix this G.D. so we don't have to always jump back into the domain
		// (causes inaccuracies in the optimization at the boundary)

//		if ( val != xval[i] ) {
//			std::cout << vars[i]->getID() << " xval " << xval[i] <<
//					" changed to " << val << std::endl;
//		}

		assert( !std::isnan( xval[i]  ) );

		if ( val != xval[i] ) res = false;

		vars[i]->assign( val );
		func.onVarAssigned( vars[i]->getID(), val );
	}

	return res;
}

} // namespace rdis
