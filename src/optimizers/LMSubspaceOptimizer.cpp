/*
 * LMSubspaceOptimizer.cpp
 *
 *  Created on: Apr 30, 2014
 *      Author: afriesen
 */

#include "optimizers/LMSubspaceOptimizer.h"

#include "util/numeric.h"

#ifdef USE_LEVMAR
#include "levmar.h"
#endif

//#include BOOSTPATH/bind.hpp>

namespace rdis {

LMSubspaceOptimizer::LMSubspaceOptimizer( OptimizableFunction & f_ )
	: super( f_ )
{
	assert( !doAscent );
}

LMSubspaceOptimizer::~LMSubspaceOptimizer() {}


Numeric LMSubspaceOptimizer::optimize( const VariablePtrVec & vars,
		const FactorPtrVec & factors, NumericVec & xval,
		Numeric & deltaFval, const bool printdbg ) {

//	// box-constrained minimization
//	extern int dlevmar_bc_der(
//	   void (*func)(double *p, double *hx, int m, int n, void *adata),
//	   void (*jacf)(double *p, double *j, int m, int n, void *adata),
//	   double *p, double *x, int m, int n, double *lb, double *ub, double *dscl,
//	   int itmax, double *opts, double *info, double *work, double *covar,
//	   void *adata);

//#define USE_LEVMAR
#ifdef USE_LEVMAR
	const Clock::time_point starttime = Clock::now();

//	initstate.assign( xval.begin(), xval.end() );
	LMSSOpt::AuxData aux( vars, factors, f, pgtemp );

	const int m = vars.size();
	const int n = std::max( factors.size(), (size_t) m );

	initeval.resize( n );
	finaleval.resize( n );

//	double lb[ m ];
//	double ub[ m ];
//
//	for ( size_t i = 0; i < vars.size(); ++i ) {
//		NumericInterval dom = vars[i]->getDomain().interval();
//		lb[i] = ( dom.lower() <= -DBL_MAX ? DBL_MAX : dom.lower() );
//		ub[i] = ( dom.upper() >= DBL_MAX ? DBL_MAX : dom.upper() );
//	}

	double * dscl = NULL;

	double opts[ LM_OPTS_SZ ];
//	double * opts = NULL;

	double info[ LM_INFO_SZ ];
	double * work = new double[ LM_BC_DER_WORKSZ( m, n ) ];

	double * covar = NULL;

	if ( printdbg ) {
		std::cout << "LM SS opt m=" << m << ", n=" << n << " (" << vars.size()
				<< " vars, " << factors.size() << " factors)" << std::endl;
	}

//	LMSSOpt::evalFunc( xval.data(), initeval.data(), m, n, &aux );
//	const Numeric ival = NumericVecOps::dot( initeval, initeval ) / 2.0;
	Numeric ferr = 0.0;
	const Numeric ival = f.evalFactors( factors, ferr );

	opts[0] = 1e-3; 	// initial \mu scale factor
	opts[1] = 1e-15;	// stopping thresh. for ||J^T e||_inf
	opts[2] = 1e-15;	// stopping thresh. for ||Dp||_2
	opts[3] = ftol; 	// stopping thresh. for ||e||_2

//	int niters = dlevmar_bc_der(
//			&LMSSOpt::evalFunc,
//			&LMSSOpt::evalJacf,
//			xval.data(), NULL, xval.size(), n, lb, ub,
//			NULL, maxiters, opts, info, work, covar, &aux );

	int niters = dlevmar_der(
			&LMSSOpt::evalFunc,
			&LMSSOpt::evalJacf,
			xval.data(), NULL, xval.size(), n,
			maxiters, opts, info, work, covar, &aux );

	// sanitize values to be within the domain of this variable (note that they
	// are sanitized in quickAssign(), so just copy them out after)
	LMSSOpt::quickAssignVals( aux, m, xval.data() );
	for ( size_t i = 0; i < vars.size(); ++i ) {
		xval[i] = vars[i]->eval();
	}

//	LMSSOpt::evalFunc( xval.data(), finaleval.data(), m, n, &aux );
//	const Numeric fval = NumericVecOps::dot( finaleval, finaleval ) / 2.0;

	Numeric fval = f.evalFactors( factors, ferr );

//	Numeric asdf = 0;
//	for ( size_t i = 0; i < factors.size(); ++i ) {
//		asdf += factors[i]->evalNoCache();
//		std::cout << i << ": f " << factors[i]->getID() << " diff " <<
//				( finaleval[i]*finaleval[i]/2.0 - factors[i]->eval() ) << " ("
//				<< factors[i]->eval() << ")" << std::endl;
//	}
//	std::cout << "asdf: " << asdf << std::endl;
//
//	std::cout << "f.eval " << f.eval() << std::endl;
//	Numeric ferr = 0.0;
//	std::cout << "f fact eval " << f.evalFactors( factors, ferr ) << std::endl;
//
//	std::cout << "xval: " << xval << std::endl;

//	double err[n];
//	dlevmar_chkjac( &LMSSOpt::evalFunc, &LMSSOpt::evalJacf, xval.data(), m, n,
//			&aux, err );
//
//	for ( size_t i = 0; i < m; ++i ) {
//		if ( err[i] < 0.5 ) {
//			std::cout << i << ": var " << vars[i]->getID() << " has jac err "
//					<< err[i] << std::endl;
//		}
//	}

	delete work;

	deltaFval = ( fval - ival);
	const Duration dur = ( Clock::now() - starttime );

	if ( printdbg ) {
		std::cout << "LM SS opt returned " << fval << " (init " << ival
				<< ", diff " << deltaFval << ") after " << niters <<
				" iterations in " << dur.count() << " seconds" << std::endl;

		std::cout << "LM SS opt info -- termination: " << info[6] <<
				", #fevals: " << info[7] << ", #jevals: " << info[8] <<
				", #linsolves" << info[9] << std::endl;
	}

//	if ( deltaFval > 0 ) {
//		fval = ival;
//		xval.assign( initstate.begin(), initstate.end() );
//		deltaFval = 0;
//		std::cout << "LM SS made no progress -- returning initial state" << std::endl;
//	}

	return fval;

#else // USE_LEVMAR

	std::cout << "ERROR: Cannot use LM subspace optimizer without levmar.h." <<
			" Rebuild with -DUSE_LEVMAR and link to liblevmar." << std::endl;

	std::exit( -1 );

	return 0;
#endif // USE_LEVMAR
}


namespace LMSSOpt {

void evalFunc( double * p, double * hx, int m, int n, void * adata_in ) {

	// void (*func)(double *p, double *hx, int m, int n, void *adata)
	//   *p is the (in) array containing the parameter values (xval)
	//   *hx (\hat(x)) is the (out) array containing \hat(x) = f(p)
	//   m is the number of parameters (length of p)
	//   n is the number of dimensions per factor (length of hx)
	//   adata_in is an auxiliary data vector used to compute \hat(x) = f(p)

	assert( p != NULL );
	assert( hx != NULL );
	assert( adata_in != NULL );

	AuxData * aux = static_cast< AuxData * >( adata_in );

	quickAssignVals( *aux, m, p );

//	Numeric ferr = 0;
//	hx[0] = aux->f.evalFactors( aux->factors, ferr );

	assert( (size_t) n >= aux->factors.size() );
	for ( int i = 0; i < n; ++i ) {
		if ( (size_t) i < aux->factors.size() ) {
			hx[i] = std::sqrt( aux->factors[i]->eval() * 2.0 );
		} else {
			hx[i] = 0;
		}
	}
}


void evalJacf( double * p, double * jac, int m, int n, void * adata_in ) {

	// void (*jacf)(double *p, double *j, int m, int n, void *adata)
	//   *p is the (in) array containing the parameter values (xval)
	//   *j is the (out) array containing the jacobian (df/dvars)
	//   m is the number of parameters (vars)
	//   n is the number of measurements (one per factor, currently)
	//   adata_in is an auxiliary data vector used to compute the jacobian

	assert( p != NULL );
	assert( jac != NULL );
	assert( adata_in != NULL );

	AuxData * aux = static_cast< AuxData * >( adata_in );

	assert( n >= aux->factors.size() );
//	assert( 2*aux->factors.size() == n );

	quickAssignVals( *aux, m, p );

//	// compute gradient
//	aux->pgt1.clear();
//	for ( Factor * f : aux->factors ) {
//		aux->pgt2.clear();
//		f->computeGradient( aux->pgt2, false );
//		productGradient( aux->pgt1, aux->pgt2 );
//	}
//
////	deriv.resize( vars.size(), 0 );
//
////	std::cout << "derivative of f() at " << x << " -- pg: ";
////	for ( auto it = pgtemp1.begin(); it != pgtemp1.end(); ++it ) {
////		std::cout << "(" << it->first << ", " << it->second << "), ";
////	}
////	std::cout << std::endl;
//
//	// copy gradient into vector
//	for ( size_t i = 0; i < aux->vars.size(); ++i ) {
//#ifdef DEBUG
//		jac[i] = aux->pgt1.at( aux->vars[i]->getID() );
//#else
//		auto it = pgtemp1.find( aux->vars[i]->getID() );
//		jac[i] = ( it == pgtemp1.end() ? 0 : it->second );
//#endif
//	}

	assert( aux->f.semiring().Product( 2, 4 ) == 6 );

	VariableID vid = -1;
	Numeric feval = 0.0;
	for ( size_t j = 0; j < n; ++j ) { // for each factor (plus padding)

		feval = 0.0;

		aux->pgt.clear();
		if ( j < aux->factors.size() ) {
			aux->factors[j]->computeGradient( aux->pgt, false );
			feval = std::sqrt( aux->factors[j]->eval() * 2.0 );
		}

		for ( size_t i = 0; i < m; ++i ) { // for each var
			vid = aux->vars[i]->getID();
			auto lb = aux->pgt.lower_bound( vid );
			if ( lb == aux->pgt.end() || lb->first < vid ) {
				jac[j*m+i] = 0;
			} else {
				jac[j*m+i] = lb->second / feval;
			}
		}
	}

}


void quickAssignVals( AuxData & aux, int m, double * p ) {
	const auto & vars = aux.vars;
	assert( vars.size() == m );
	for ( size_t i = 0; i < vars.size(); ++i ) {
		Numeric val = p[i];
//		Numeric val = vars[i]->getDomain().closestVal( p[i] );
		// TODO: Should we ensure we're within the domain here?


//		if ( val != p[i] ) {
//			std::cout << vars[i]->v->getID() << " xval " << p[i] << " changed to " << val << std::endl;
//		}
//		Numeric val = p[i];

		vars[i]->assign( val );
		aux.f.onVarAssigned( vars[i]->getID(), val );
	}
}

} // namespace LMSSOpt

} // namespace rdis
