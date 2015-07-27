/*
 * BCDOptimizer.cpp
 *
 *  Created on: Jan 26, 2014
 *      Author: afriesen
 */

#include "optimizers/BCDOptimizer.h"

#include BOOSTPATH/chrono.hpp>

namespace rdis {

BCDOptimizer::BCDOptimizer( OptimizableFunction & f_, SubspaceOptimizer & ssopt_,
		size_t blocksize_, Numeric epsilon_ )
	: fnc( f_ )
	, ssopt( ssopt_ )
	, blocksize( blocksize_ )
	, epsilon( epsilon_ )
	, timedOut( false )
	, fopt( std::numeric_limits< Numeric >::max() )
	, numIterations( 0 )
	, rng( 4582237 )
{
}

BCDOptimizer::~BCDOptimizer() {}


Numeric BCDOptimizer::optimize( long long int maxitersCG,
		long long int maxitersTotal, const State * xinitial,
		bool printdbg, Clock::time_point timeout, size_t numRandomRestarts,
		Clock::time_point actualStartTime ) {

	Numeric result( 0 ), prevResult( std::numeric_limits< Numeric >::max() );
	const Clock::time_point starttime = ( actualStartTime < Clock::now() ) ?
				actualStartTime : Clock::now();

	timedOut = false;
	xopt.clear();
	xopt.resize( fnc.getNumVars(), 0 );

	VariablePtrVec & vars( fnc.getVariables() );
	State x( fnc.getNumVars(), 0 );

	NumericVec gradient;

	for ( VariableID vid = 0; vid < fnc.getNumVars(); ++vid ) {
		x[vid] = ( xinitial == NULL ? vars[vid]->getDomain().median() :
				(*xinitial)[vid] );
	}

	size_t nblocks = std::ceil( (Numeric) fnc.getNumVars() ) /
			( (Numeric) blocksize );
	blocks.resize( nblocks + 2 );
	fblocks.resize( nblocks + 2 );

	// create the blocks

	size_t bid = 0;
	VariableID vbl, vbu;
	for ( VariableCount vbid = 0; vbid < fnc.getNumBlocks(); ++vbid ) {
		fnc.getBlockRangeByBlkId( vbid, vbl, vbu );
		if ( vbl < 0 || vbu < 0 ) continue;

		const VariableCount vbsz = vbu - vbl + 1;
		assert( (VariableCount) blocksize >= vbsz ); // block size can't be too small

		if ( blocks[bid].size() + vbsz > blocksize+1 ) { // +1 is to allow squeezing in one extra if necessary
			++bid;
			assert( blocks.size() > bid );
			assert( blocks[bid].size() + vbsz <= blocksize );
		}

		for ( VariableID vid = vbl; vid <= vbu; ++ vid ) {
			blocks[bid].push_back( vars[vid] );
			fblocks[bid].insert( fblocks[bid].end(),
					vars[vid]->getFactors().begin(),
					vars[vid]->getFactors().end() );

//			std::cout << "putting var " << vid << "(" << vbid << ") in block " << bid << std::endl;
		}
	}

	// remove any empty blocks (caused by e.g. floating point noise)
	while ( blocks.back().empty() ) {
		blocks.pop_back();
		fblocks.pop_back();
	}

	nblocks = blocks.size();

	// remove all duplicates from the factor lists
	for ( auto & fb : fblocks ) {
		std::sort( fb.begin(), fb.end() );
		auto it = std::unique( fb.begin(), fb.end() );
		fb.erase( it, fb.end() );
	}

	std::cout << "there are " << nblocks << " blocks (" << blocksize << ")" << std::endl;
	for ( int i= 0; i < (int) blocks.size(); ++i ) {
		std::cout << "block " << i << " has " << blocks[i].size() << " vars and "
				<< fblocks[i].size() << " factors" << std::endl;
	}

	VariableCount iteration = 0;

	State xtemp;
	Numeric mineval = std::numeric_limits< Numeric >::max();

	for ( size_t rrID = 0; rrID <= numRandomRestarts; ++rrID ) {

		// generate a random starting state
		if ( rrID > 0 ) {
			const NumericInterval dom = vars[0]->getDomain().getSamplingInterval();
			BOOSTNS::random::uniform_real_distribution<> unif( dom.lower(), dom.upper() );

			for ( VariableID vid = 0; vid < (VariableID) vars.size(); ++vid ) {
				x[vid] = unif( rng );
			}

			std::cout << "BCD converged, so randomly restarting" << std::endl;
		}

		assignVars( x );

		// compute the initial log likelihood
		prevResult = result = fnc.eval();

		for ( ; ; ) {

			if ( printdbg /*&& iteration == 0*/ ) {
				std::cout << "iteration " << iteration << ", eval: " << result <<
						", elapsed: " << (Duration) ( Clock::now() - starttime ) << std::endl;
				std::cout << "\tstate: " << x << std::endl;
			}

			// do GD on each block
			for ( size_t bid = 0; bid < nblocks; ++bid ) {
				size_t vidx = 0;
				xtemp.resize( blocks[bid].size() );

				// copy initial state in
				for ( Variable * v : blocks[bid] ) {
					xtemp[vidx++] = x[v->getID()];
				}

				Numeric deltafval = 0.0;
				ssopt.optimize( blocks[bid], fblocks[bid], xtemp, deltafval, printdbg );

				vidx = 0;

				// copy final state out
				for ( Variable * v : blocks[bid] ) {
					x[v->getID()] = xtemp[vidx++];
				}

				std::cout << "NEW FULL EVAL: " << fnc.eval() << " after " <<
						( (Duration) ( Clock::now() - starttime ) ).count() <<
						" seconds" << std::endl;

				if ( Clock::now() > timeout ) {
					timedOut = true;
					break;
				}
			}

			// assign the new values
			assignVars( x );

			// compute the new log likelihood
			prevResult = result;
			result = fnc.eval();

			std::cout << "new eval: " << result << std::endl;

			if ( std::isnan( result ) ) {
				if ( printdbg ) std::cerr << "BCD returned NAN" << std::endl;
				result = prevResult;
				break;
			}

			if ( result < mineval ) {
				mineval = result;
				xopt = x;
			}

			// if we've converged, get the final state and unassign all vars
			if ( ++iteration >= maxitersTotal && maxitersTotal >= 0 ) break;
			if ( -( result - prevResult ) < epsilon ) break;

			// check for timeout
			if ( Clock::now() > timeout ) {
				std::cout << "BCD opt timed out" << std::endl;
				break;
			}
		}

		if ( Clock::now() > timeout ) break;
	}

	if ( printdbg ) {
		std::cout << "iteration " << iteration << ": eval: " << result << std::endl;
	}

	for ( VariableID vid = 0; vid < fnc.getNumVars(); ++vid ) {
		vars[vid]->unassign();
		fnc.onVarUnassigned( vid );
	}


	fopt = mineval;
	runtime = ( Clock::now() - starttime );
	numIterations = iteration;

	if ( printdbg ) {
		std::cout << "**********" << std::endl;
		std::cout << "Block coordinate descent completed in " << iteration <<
				" iterations and " << runtime.count() << " seconds" << std::endl;
		std::cout << "\tresult: " << fopt << std::endl;
		std::cout << "\tstate: " << xopt << std::endl;
		std::cout << "**********" << std::endl;
	}

	return fopt;
}


void BCDOptimizer::assignVars( const State & x ) {
	VariablePtrVec & vars = fnc.getVariables();
	for ( VariableID vid = 0; vid < fnc.getNumVars(); ++vid ) {
		Numeric val = vars[vid]->getDomain().closestVal( x[vid] );
		vars[vid]->assign( val );
		fnc.onVarAssigned( vid, val );
	}
}

} // namespace rdis
