/*
 * EMOptimizer.cpp
 *
 *  Created on: Sep 9, 2013
 *      Author: afriesen
 */

#include "optimizers/EMOptimizer.h"

#include BOOSTPATH/math/distributions/normal.hpp>


namespace rdis {

EMOptimizer::EMOptimizer( GMMFunction & fgmm_, Numeric epsilon_ )
	: fgmm( fgmm_ )
	, epsilon( epsilon_ )
	, fopt( std::numeric_limits< Numeric >::max() )
	, runtime( 0 )
	, numIterations( 0 )
{
}

EMOptimizer::~EMOptimizer() {}



Numeric EMOptimizer::optimize( long long int maxiters, const State * xinitial,
		bool printdbg, Clock::time_point /*timeout*/ ) {

	Numeric initresult = 0;
	Numeric result( 0 ), prevResult( std::numeric_limits< Numeric >::max() );
	const Clock::time_point starttime( Clock::now() );

	xopt.clear();
	xopt.resize( fgmm.getNumVars(), 0 );

	// TODO: assume randomly initialized?? or do it here? or use xinitial?

	VariablePtrVec & vars( fgmm.getVariables() );

	// assign all variables their initial values
	for ( VariableID vid = 0; vid < fgmm.getNumVars(); ++vid ) {
		vars[vid]->assign( xinitial == NULL ?
						   vars[vid]->getDomain().median() :
				(*xinitial)[vid] );
	}

	// compute the initial log likelihood
	initresult = prevResult = result = computeLL();

	long long int iteration = 0;
	bool failed = false;

	for ( ; ; ) {

//		if ( printdbg && iteration == 0 ) {
//			std::clog << "iteration " << iteration << ": neg. log likelihood: " <<
//					-result << std::endl;
//		}

		// expectation step: compute Q(\th|\th^{p}) = E[ log( f(x|\th') ) | y, \th ]
		EStep();

		// maximization step: compute \th^{p+1} = argmax_{\th} Q(\th|\th^{p})
		if ( !MStep() ) {
			prevResult = result;
			result = computeLL();
			if ( printdbg ) {
				std::cerr << "iteration " << iteration <<
						": EM failed in M step at " << -result << " (" <<
						-prevResult << ")" << std::endl;
			}
			failed = true;
			break;
		}

//		fgmm.EMStep( vars, fgmm.getFactors(), false );

		// compute the new log likelihood
		prevResult = result;
		result = computeLL();

		// if we've converged, get the final state and unassign all vars
		if ( ++iteration >= maxiters && maxiters >= 0 ) break;
		if ( ( result - prevResult ) < epsilon ) break;

//		break;
	}

//	if ( printdbg ) {
//		std::clog << "iteration " << iteration << ": neg. log likelihood: "
//				<< -result << std::endl;
//	}

	for ( VariableID vid = 0; vid < fgmm.getNumVars(); ++vid ) {
		xopt[vid] = vars[vid]->eval();
		vars[vid]->unassign();
	}

	fopt = result;
	runtime = ( Clock::now() - starttime );
	numIterations = iteration;

	if ( printdbg ) {
//		std::cout << "**********" << std::endl;
//		std::cout << "EM optimization completed in " << iteration <<
//				" iterations and " << runtime.count() << " seconds" << std::endl;
//		std::cout << "\tneg. log likelihood: " << -result << std::endl;
//		std::cout << "\tstate: "; printVector( xopt ); std::cout << std::endl;
//		std::cout << "**********" << std::endl;
//	} else {
		std::cout << "EM opt. (iters " << iteration << "): " << -result
				<< " (init " << -initresult << ", diff "
				<< ( initresult - result ) << ") in " << runtime.count()
				<< " seconds" << std::endl;
//		std::cout << "iters? " << ( iteration >= maxiters && maxiters >= 0 )
//				<< ", convergence? " << ( ( result - prevResult ) < epsilon )
//				<< std::endl;
	}

	return result;
}


void EMOptimizer::EStep() {
	// p_{i,k} = a_k^(p) * f(x_i | \mu_k^(p), \sigma_k^(p) ) / ...
	// 				\sum_{k=1}^{K} a_k^(p) f(x_i | \mu_k^(p), \sigma_k^(p) )

	const size_t N = fgmm.getNumObservations();
	const size_t K = fgmm.getNumClusters();

	assert( fgmm.getFactors().size() == N );
	posteriorProbs.resize( N, std::vector< Numeric >( K ) );

	// TODO: speed-up this computation? make it consistent with the GMMFactor
	// computation?

	for ( size_t i = 0; i < N; ++i ) {

		posteriorProbs[i].clear();
		posteriorProbs[i].resize( K, 0 );

		Numeric Z( 0.0 );

		// compute the unnormalized posterior probabilities
		for ( size_t k = 0; k < K; ++k ) {
			posteriorProbs[i][k] = feval( i, k );
			assert( !std::isnan( posteriorProbs[i][k] ) );
			Z += posteriorProbs[i][k];
		}

		// normalize all of the posterior probabilities
		for ( auto & p_ik : posteriorProbs[i] ) {
			if ( Z == 0 ) {
				assert( p_ik == 0 );
			} else if ( std::isinf( Z ) ) {
				p_ik = 0;
			} else {
				p_ik /= Z;
			}
			assert( !std::isnan( p_ik ) );
		}
	}

//	std::cout << "p_ik (o): " << std::endl;
//	for ( size_t k = 0; k < K; ++k ) {
//		std::cout << "cluster " << k << ": ";
//		for ( size_t i = 0; i < N; ++i ) {
//			std::cout << posteriorProbs[i][k] << ", ";
//		}
//		std::cout << std::endl;
//	}
}



bool EMOptimizer::MStep() {
	if ( !maximizeWeights() ) return false;
	if ( !maximizeMu() ) return false;
	if ( !maximizeCov() ) return false;
	return true;
}


Numeric EMOptimizer::computeLL() {
	const size_t N = fgmm.getNumObservations();
	const size_t K = fgmm.getNumClusters();

	Numeric loglik( 0 );

	for ( size_t i = 0; i < N; ++i ) {
		Numeric lli = 0;
		for ( size_t k = 0; k < K; ++k ) {
			lli += feval( i, k );
		}
		loglik += std::log( lli );
	}

	return loglik;
}


Numeric EMOptimizer::feval( size_t i, size_t k ) {
//	const Numeric c( -0.918938533204672741780329736 ); // == log( 1/sqrt(2*pi))
	const Numeric c = -0.91893853320467274178032973640561763986139747363778; // == ln(1/sqrt(2*pi))
//	const Numeric c( std::log( 1.0 /
//			BOOSTNS::math::constants::root_two_pi< Numeric >() ) );
	Numeric logteval( 0.0 );

	const Numeric D( fgmm.getNumDims() );

	const VariableID meanID( fgmm.getFirstMu( k ) );
	const VariableID covID( fgmm.getFirstCov( k ) );
	const VariableID wtID( fgmm.getWeight( k ) );

	const VariablePtrVec & vars = fgmm.getVariables();
	const NumericVecVec & obs = fgmm.getObservations();

	for ( size_t d = 0; d < D; ++d ) {
		Numeric mu = vars[meanID+d]->eval();
		Numeric sig2 = vars[covID+d]->eval();

		Numeric logp = c - std::log( sig2 ) / 2 -
				( obs[i][d] - mu ) * ( obs[i][d] - mu ) / ( 2.0 * sig2 );

		assert( !std::isnan( logp ) );
//		assert( std::exp( p ) <= 1.0 );
		logteval += logp;

//#ifdef DEBUG
//		// TODO: remove this check
//		BOOSTNS::math::normal n( mu, std::sqrt( sig2 ) );
//		Numeric p2 = BOOSTNS::math::pdf( n, obs[i][d] );
//		if ( !approxeq( std::exp( logp ), p2, 1e-12 ) ) {
//			std::cout << "diff: " << ( std::exp( logp ) - p2 ) << std::endl;
//			assert( approxeq( std::exp( logp ), p2, 1e-8 ) );
//		}
//#endif
	}

	Numeric eval = std::exp( logteval );
	assert( !std::isnan( eval ) );
//	assert( eval <= 0.0 );

	return eval * vars[wtID]->eval();
}


bool EMOptimizer::maximizeWeights() {
	const size_t N = fgmm.getNumObservations();
	const size_t K = fgmm.getNumClusters();

	for ( size_t k = 0; k < K; ++k ) {
		Numeric wnew = 0;
		for ( size_t i = 0; i < N; ++i ) {
			wnew += posteriorProbs[i][k];
		}
		wnew /= N;

#ifndef EM_UPDATEWEIGHTS
		wnew = 1.0 / ( (Numeric) K );
#endif

		if ( wnew <= 0.0 ) return false;

		fgmm.getVariables()[fgmm.getWeight( k )]->assign( wnew );

//		std::cout << "setting weight var " << fgmm.getWeight(k) << " to " << wnew << std::endl;
	}

	return true;
}


bool EMOptimizer::maximizeMu() {
	const size_t N = fgmm.getNumObservations();
	const size_t K = fgmm.getNumClusters();
	const size_t D = fgmm.getNumDims();

	const NumericVecVec & obs = fgmm.getObservations();

	for ( size_t k = 0; k < K; ++k ) {
		for ( size_t d = 0; d < D; ++d ) {
			Numeric xpik = 0;
			Numeric pik = 0;
			for ( size_t i = 0; i < N; ++i ) {
				xpik += obs[i][d] * posteriorProbs[i][k];
				pik += posteriorProbs[i][k];
			}

			if ( pik == 0.0 ) return false;

			Numeric munew = xpik / pik;
			fgmm.getVariables()[fgmm.getFirstMu( k ) + d]->assign( munew );

//			std::cout << "setting mu var " << fgmm.getFirstMu( k ) + d << " to " << munew << std::endl;
		}
	}

	return true;
}


bool EMOptimizer::maximizeCov() {
	const size_t N = fgmm.getNumObservations();
	const size_t K = fgmm.getNumClusters();
	const size_t D = fgmm.getNumDims();

	VariablePtrVec & vars = fgmm.getVariables();
	const NumericVecVec & obs = fgmm.getObservations();

	// Note that the following computation uses the just-computed value of mu,
	// so this function must be called after maximizeMu() in MStep().

	for ( size_t k = 0; k < K; ++k ) {
		for ( size_t d = 0; d < D; ++d ) {
			Numeric pikxmu = 0;
			Numeric pik = 0;
			for ( size_t i = 0; i < N; ++i ) {
				Variable * vmu = vars[fgmm.getFirstMu( k ) + d];
				pikxmu += posteriorProbs[i][k] *
						( obs[i][d] - vmu->eval() ) * ( obs[i][d] - vmu->eval() );
				pik += posteriorProbs[i][k];
			}
			Numeric covnew = 0.0;
			if ( pik == 0.0 ) {
				assert( pikxmu == 0.0 );
			} else {
				covnew = pikxmu / pik;
			}

			if ( covnew <= 0.0 ) return false;
			assert( covnew > 0.0 );
			vars[fgmm.getFirstCov( k ) + d]->assign( covnew );

//			std::cout << "setting cov var " << fgmm.getFirstCov( k ) + d << " to " << covnew << std::endl;
		}
	}

	return true;
}

} // namespace rdis
