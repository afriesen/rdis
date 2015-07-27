/*
 * optimize_sinusoid.cpp
 *
 *  Created on: Jan 26, 2014
 *      Author: afriesen
 */


#include "common.h"
#include "OptimizableFunction.h"
#include "OptimizableFunctionGenerator.h"

#include "optimizers/CGDSubspaceOptimizer.h"

#include "RDISOptimizer.h"

#include "optimizers/BCDOptimizer.h"

#include "util/utility.h"
#include "util/signal_handler.h"

#include BOOSTPATH/foreach.hpp>
#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/filesystem.hpp>

#include BOOSTPATH/assign/list_of.hpp> // for 'list_of()'
#include BOOSTPATH/assign/std/vector.hpp> // for 'operator+=()'

#include BOOSTPATH/math/distributions/normal.hpp>


namespace fs = BOOSTNS::filesystem;
namespace bc = BOOSTNS::chrono;
namespace psi = rdis::instrumentation;
namespace po = BOOSTNS::program_options;

using namespace rdis;



void optimize_sinusoid( BOOSTNS::random::mt19937 & rng,
						const po::variables_map & options );

bool doBCDOpt( OptimizableFunction & func, SubspaceOptimizer & ssopt,
		const po::variables_map & options, const State & xinit,
		Clock::time_point timeout, Clock::time_point starttime,
		State & xopt, Numeric & fopt );

Numeric eval( OptimizableFunction & f, const State & x );
void addOptions( po::options_description & desc );


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main( int argc, const char * const argv[] )
{
	BOOSTNS::random::mt19937 rng( 834725927 );
	const Clock::time_point start = Clock::now();

	std::cout << "called with: " << std::endl << "\t";
	for ( int i = 0; i < argc; ++i ) {
		std::cout << argv[i] << " ";
	}
	std::cout << std::endl;

	po::options_description desc( "Sinusoid optimization options" );

	try {
		addOptions( desc );

		po::variables_map options;
		po::store( po::parse_command_line( argc, argv, desc ), options );

		po::notify( options );

		if ( options.count( "help" ) ) {
			std::cout << desc << std::endl;
		} else {
			optimize_sinusoid( rng, options );
		}

	} catch ( po::error & e) {
		  std::cerr << "\nCommand line error: " << e.what() << std::endl << std::endl;
		  std::cerr << desc << std::endl;
	} catch ( std::exception & e ) {
		std::cerr << "\nException occurred in main (p.f. opt): \n\t" << e.what()
				<< ", application will now exit. Goodbye!" << std::endl;
//	} catch ( utility::excn::EXCN_Base & e ) {
//		std::cerr << "\nCaught rosetta exception: " << e.msg() <<
//				", application will now exit. Goodbye!" << std::endl;
	}

	const Duration elapsed = Clock::now() - start;
	std::cout << "total elapsed time: " << elapsed.count() << " seconds" << std::endl;
}



void optimize_sinusoid( BOOSTNS::random::mt19937 & rng,
						const po::variables_map & options ) {

	Numeric regionwidth = 1.0;
	Numeric eps = 0.1, lipsch = 20;
	RDISOptimizer * popt = NULL;

	const string domain =
			BOOSTNS::str( BOOSTNS::format( "%1%:%2%" ) %
					( -regionwidth/2.0 ) % ( regionwidth/2.0 ) );
//					( regionwidth / 4.0 ) % ( 5.0*regionwidth/4.0 ) );

	eps = options.count( "epsilon" ) ? options[ "epsilon" ].as< Numeric >() : eps;
	lipsch = options.count( "lipschitz" ) ? options[ "lipschitz" ].as< Numeric >() : lipsch;

	if ( options.count( "timeAsSeed" ) && options[ "timeAsSeed" ].as< bool >() ) {
		rng.seed( std::time(NULL) + getpid() );
		if ( options[ "sampleOffset"].as< size_t >() > 0 ) {
			std::cout << "WARNING: sampleOffset and timeAsSeed should not both be set" << std::endl;
		}
	}

	const bool isBCD = ( options.count( "opttype" ) &&
			options[ "opttype" ].as< string >() == "bcd" );

	size_t h = options[ "sinHeight" ].as< size_t >();
	size_t bf = options[ "sinBF" ].as< size_t >();
	size_t arity = options[ "sinArity" ].as< size_t >();
	bool oddArity = options[ "sinOddArity" ].as< bool >();

	OptimizableFunctionSP funcsp =
			OptimizableFunctionGenerator::makeHighDimSinusoid( h, bf, arity, oddArity );

	OptimizableFunction & fsinusoid = *funcsp;

	std::cout << fsinusoid << std::endl;

	CGDSubspaceOptimizer ssopt( fsinusoid );
	ssopt.setParameters( options );

	// create the optimizer
	popt = new RDISOptimizer( fsinusoid, ssopt, options );

	RDISOptimizer & opt = *popt;
	setSignalHandler( &opt );

	opt.setApproxFactors( options[ "approxfactors" ].as< bool >() );
	opt.setUseCache( false );
//	opt.setUseExactCaching( false );

	const Numeric MaxVal = std::numeric_limits< Numeric >::max();

	size_t sampleOffset = options[ "sampleOffset" ].as< size_t >();
	size_t nsamples = options[ "nsamples" ].as< size_t >() + sampleOffset;

	const NumericInterval dom = fsinusoid.getVariables()[0]->getDomain().interval();

	BOOSTNS::random::uniform_real_distribution<> unif( dom.lower(), dom.upper() );
	std::vector< State > initialStates( nsamples );

	// create all the random initial states up front (so they're consistent
	// across runs if the seed doesn't change)
	for ( size_t i = 0; i < nsamples; ++i ) {
		for ( VariableID vid = 0; vid < fsinusoid.getNumVars(); ++vid ) {
			initialStates[i].push_back( unif( rng ) );
		}
	}

	const Clock::time_point start = Clock::now();
	Clock::time_point timeout = Clock::time_point::max();

	if ( options.count( "timeout" ) ) {
		size_t ms = 1000.0 * options[ "timeout" ].as< Numeric >();
		timeout = Clock::now() + bc::milliseconds( ms );
	}

	Numeric bestfopt = NumericMAX;
	State bestx;

	size_t i = sampleOffset;
	for ( ; i < nsamples; ) {

		Numeric finit = MaxVal;
		const State & xinit = initialStates[i];

		std::cout << "[" << i << "] optimization starting -- bounds: " <<
				fsinusoid.computeBounds() << std::endl;
		std::cout << "eps: " << eps << ", lipschitz: " << lipsch <<
				", step: " << eps / lipsch <<
				", finit: " << finit << std::endl << "xinit: " << xinit << std::endl;

		std::cout << "initial eval: " << eval( fsinusoid, xinit ) << std::endl;

		++i;

		if ( options.count( "useDA" ) ) {
			throw "deterministic annealing not yet supported";
		} else if ( options.count( "useCP" ) ) {
			throw "cutting plane not yet supported";
		}

		State xopt;
		Numeric result = MaxVal;
		bool timedOut = false;

		if ( !isBCD ) {
			result = opt.optimize( &xinit, finit, timeout, true, start );
			xopt = opt.getOptState();
			opt.printStats();
			timedOut = ( opt.wasStopped() || opt.timedOut() );
		} else {
			timedOut = doBCDOpt( fsinusoid, ssopt, options, xinit, timeout, start,
					xopt, result );
		}

		const Duration elapsed = Clock::now() - start;

		assert( !xopt.empty() );
		std::cout << "sinusoid opt. result: " << result << std::endl <<
					"\tat " << xopt << std::endl;
		Numeric actual = eval( fsinusoid, xopt );
		std::cout << "actual eval: " << actual << " after " << elapsed.count()
				<< " seconds " << std::endl;

		if ( actual < bestfopt ) {
			bestfopt = actual;
			bestx = xopt;
		}

		if ( timedOut ) {
			std::cout << "timeout detected..." << std::endl;
			break; break;
		}
	}
}


bool doBCDOpt( OptimizableFunction & func, SubspaceOptimizer & ssopt,
		const po::variables_map & options, const State & xinit,
		Clock::time_point timeout, Clock::time_point starttime,
		State & xopt, Numeric & fopt ) {

	size_t blksize = std::ceil( (Numeric) func.getNumVars() / 5.0 );
	if ( options.count( "AVblksz" ) ) {
		blksize = options["AVblksz"].as< size_t >();
	} else if ( options.count( "AVblkpct" ) ) {
		blksize = round( options["AVblkpct"].as< Numeric >() *
				(Numeric) func.getNumVars() );
	}

	Numeric tol = 1e-6;
	if ( options.count( "steptol" ) ) tol = options["steptol"].as< Numeric >();

	size_t cgmaxit = 100;
	if ( options.count( "SSmaxit" ) ) cgmaxit = options["SSmaxit"].as< VariableCount >();

	size_t nRR = 0;
	if ( options.count( "nRRatTop" ) ) nRR = options["nRRatTop"].as< size_t >();

	std::cout << "BCD OPTIMIZATION" << std::endl;

	BCDOptimizer bcdopt( func, ssopt, blksize, tol );

	fopt = bcdopt.optimize( cgmaxit, 100000, &xinit, true, timeout, nRR,
			starttime );

	xopt = bcdopt.getOptState();
	std::cout << "BCD opt eval " << fopt << " in " << bcdopt.getRuntime() << std::endl;

	return bcdopt.getTimedOut();
}


Numeric eval( OptimizableFunction & f, const State & x ) {
	return f.eval( x );
}


// add all options to boost program_options
void addOptions( po::options_description & desc ) {

	desc.add_options()
//			( "file,f", po::value< std::vector< string > >(),
//					"input file (PDB format) containing protein(s) to optimize (can be multiple)" )
			( "opttype", po::value< string >()->default_value( "rdis" ),
					"type of optimization to perform (can be: rdis, bcd" )
			( "nsamples", po::value< size_t >()->default_value( 1 ),
					"number of times to optimize each specified PDB file" )
			( "sampleOffset", po::value< size_t >()->default_value( 0 ),
					"sample number to start at (0-indexed)" )
//			( "perruntimeout", po::value< Numeric >(),
//					"per-run optimization timeout, in seconds (applies to det. annealing and cutting plane)" )
			( "sinHeight", po::value< size_t >()->default_value( 4 ),
					"the height of the tree-like sinusoid to generate" )
			( "sinBF", po::value< size_t >()->default_value( 3 ),
					"the branching factor of the tree-like sinusoid to generate" )
			( "sinArity", po::value< size_t >()->default_value( 3 ),
					"the maximum arity of the factors in the tree-like sinusoid to generate" )
			( "sinOddArity", po::value< bool >()->default_value( false ),
					"true to allow odd arity factors in the generated tree-like sinusoid" )
					;

		RDISOptimizer::addOptions( desc );
}
